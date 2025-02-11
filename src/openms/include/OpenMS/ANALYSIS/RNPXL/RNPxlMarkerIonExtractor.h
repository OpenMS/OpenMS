// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>

#include <map>
#include <vector>

namespace OpenMS
{
  class String;

  struct OPENMS_DLLAPI RNPxlMarkerIonExtractor
  {
    /// name to mass-intensity pair
    typedef std::map<String, std::vector<std::pair<double, double> > > MarkerIonsType;
  
    /// extract an annotate RNA marker ions
    static MarkerIonsType extractMarkerIons(const PeakSpectrum& s, const double marker_tolerance);
  };
}


