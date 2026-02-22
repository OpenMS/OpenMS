// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Sachsenberg $
// $Authors: Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/DigestionEnzymeDataProvider.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeRNA.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/config.h>

namespace OpenMS
{
  /**
    @ingroup Chemistry

    @brief Data provider that loads digestion enzymes from an XML file

    Parses enzyme definitions from a ParamXML-formatted file (e.g. Enzymes.xml,
    Enzymes_RNA.xml). The file is located via File::find().

    Template parameter @p EnzymeType should be a subclass of DigestionEnzyme
    (e.g. DigestionEnzymeProtein, DigestionEnzymeRNA).
  */
  template <typename EnzymeType>
  class EnzymeXMLDataProvider : public DigestionEnzymeDataProvider<EnzymeType>
  {
  public:
    /**
     * @brief Construct a provider for the given enzyme XML file
     *
     * @param[in] filename Relative filename to locate via File::find() (e.g. "CHEMISTRY/Enzymes.xml")
     * @param[in] optional If true, return empty vector when file is not found instead of throwing
     */
    explicit EnzymeXMLDataProvider(const String& filename, bool optional = false);

    /// @brief Parses enzyme definitions from the configured XML file.
    std::vector<std::unique_ptr<EnzymeType>> loadEnzymes() override;

  private:
    String filename_;
    bool optional_;
  };

  // Explicit instantiation declarations — definitions are in EnzymeXMLDataProvider.cpp
  extern template class OPENMS_DLLAPI EnzymeXMLDataProvider<DigestionEnzymeProtein>;
  extern template class OPENMS_DLLAPI EnzymeXMLDataProvider<DigestionEnzymeRNA>;

} // namespace OpenMS
