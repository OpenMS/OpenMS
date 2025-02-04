// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once


#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>
#include <set>
#include <map>

namespace OpenMS
{

  class MassExplainer;
  class FeatureMap;
  class ChargePair;

  class OPENMS_DLLAPI ILPDCWrapper
  {

public:
    typedef std::vector<ChargePair> PairsType;
    typedef PairsType::size_type PairsIndex;

    ///Constructor
    ILPDCWrapper();

    ///Destructor
    virtual ~ILPDCWrapper();

    /// Compute optimal solution and return value of objective function
    /// If the input feature map is empty, a warning is issued and -1 is returned.
    /// @return value of objective function
    /// and @p pairs will have all realized edges set to "active"
    double compute(const FeatureMap& fm, PairsType& pairs, Size verbose_level) const;

private:

    /// slicing the problem into subproblems
    double computeSlice_(const FeatureMap& fm,
                         PairsType& pairs,
                         const PairsIndex margin_left,
                         const PairsIndex margin_right,
                         const Size verbose_level) const;

    /// slicing the problem into subproblems
    double computeSliceOld_(const FeatureMap& fm,
                            PairsType& pairs,
                            const PairsIndex margin_left,
                            const PairsIndex margin_right,
                            const Size verbose_level) const;

    /// calculate a score for the i_th edge
    double getLogScore_(const PairsType::value_type& pair, const FeatureMap& fm) const;

    typedef std::map<String, std::set<Size> > FeatureType_;

    // add another charge annotation variant for a feature
    void updateFeatureVariant_(FeatureType_& f_set, const String& rota_l, const Size& v) const;



  }; // !class

} // !namespace

