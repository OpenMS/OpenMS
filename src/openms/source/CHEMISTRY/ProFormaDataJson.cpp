// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ProFormaDataJson.h>

namespace OpenMS
{

  String toJSON(const Peptidoform& pf)
  {
    nlohmann::json j = pf;
    return String(j.dump());
  }

  Peptidoform peptidoformFromJSON(const String& json_str)
  {
    nlohmann::json j = nlohmann::json::parse(static_cast<std::string>(json_str));
    return j.get<Peptidoform>();
  }

  String toJSON(const PeptidoformIon& pfi)
  {
    nlohmann::json j = pfi;
    return String(j.dump());
  }

  PeptidoformIon peptidoformIonFromJSON(const String& json_str)
  {
    nlohmann::json j = nlohmann::json::parse(static_cast<std::string>(json_str));
    return j.get<PeptidoformIon>();
  }

} // namespace OpenMS
