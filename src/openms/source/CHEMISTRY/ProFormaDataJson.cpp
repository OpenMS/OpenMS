// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ProFormaDataJson.h>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{

  String toJSON(const Peptidoform& pf)
  {
    nlohmann::json j = pf;
    return String(j.dump());
  }

  Peptidoform peptidoformFromJSON(const String& json_str)
  {
    try
    {
      nlohmann::json j = nlohmann::json::parse(static_cast<std::string>(json_str));
      return j.get<Peptidoform>();
    }
    catch (const nlohmann::json::exception& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        json_str, String("JSON parsing failed: ") + e.what());
    }
    catch (const std::exception& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        json_str, String("JSON deserialization failed: ") + e.what());
    }
  }

  String toJSON(const PeptidoformIon& pfi)
  {
    nlohmann::json j = pfi;
    return String(j.dump());
  }

  PeptidoformIon peptidoformIonFromJSON(const String& json_str)
  {
    try
    {
      nlohmann::json j = nlohmann::json::parse(static_cast<std::string>(json_str));
      return j.get<PeptidoformIon>();
    }
    catch (const nlohmann::json::exception& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        json_str, String("JSON parsing failed: ") + e.what());
    }
    catch (const std::exception& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        json_str, String("JSON deserialization failed: ") + e.what());
    }
  }

} // namespace OpenMS
