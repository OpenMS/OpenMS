// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/AbsoluteQuantitationStandardsFile.h>

namespace OpenMS
{
  void AbsoluteQuantitationStandardsFile::load(
    const String& filename,
    std::vector<AbsoluteQuantitationStandards::runConcentration>& run_concentrations
  ) const
  {
    CsvFile csv(filename);
    StringList sl;
    std::map<String, Size> headers;
    if (csv.rowCount() > 0) // avoid accessing a row in an empty file
    {
      csv.getRow(0, sl);
    }
    for (Size i = 0; i < sl.size(); ++i)
    {
      headers[sl[i]] = i; // for each header found, assign an index value to it
    }
    run_concentrations.clear();
    for (Size i = 1; i < csv.rowCount(); ++i)
    {
      csv.getRow(i, sl);
      run_concentrations.push_back(extractRunFromLine_(sl, headers));
    }
  }

  AbsoluteQuantitationStandards::runConcentration AbsoluteQuantitationStandardsFile::extractRunFromLine_(
    const StringList& line,
    const std::map<String, Size>& headers
  ) const
  {
    AbsoluteQuantitationStandards::runConcentration rc;
    std::map<String, Size>::const_iterator it;
    it = headers.find("sample_name");
    rc.sample_name = it != headers.end() ? line[it->second] : "";
    it = headers.find("component_name");
    rc.component_name = it != headers.end() ? line[it->second] : "";
    it = headers.find("IS_component_name");
    rc.IS_component_name = it != headers.end() ? line[it->second] : "";
    it = headers.find("actual_concentration");
    rc.actual_concentration = it != headers.end() ? line[it->second].toDouble() : 0.0;
    it = headers.find("IS_actual_concentration");
    rc.IS_actual_concentration = it != headers.end() ? line[it->second].toDouble() : 0.0;
    it = headers.find("concentration_units");
    rc.concentration_units = it != headers.end() ? line[it->second] : "";
    it = headers.find("dilution_factor");
    rc.dilution_factor = it != headers.end() ? line[it->second].toDouble() : 1.0;
    return rc;
  }
}
