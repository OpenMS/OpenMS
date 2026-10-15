// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

#include <vector>
#include <map>
#include <fstream>

namespace OpenMS
{
  namespace Internal
  {
    /**
        @brief Handler that is used for parsing PTMXML data
    */
    class OPENMS_DLLAPI PTMXMLHandler :
      public XMLHandler
    {
public:
      /// Constructor for loading
      PTMXMLHandler(std::map<std::string, std::pair<std::string, std::string> > & ptm_informations, const std::string & filename);

      /// Destructor
      ~PTMXMLHandler() override;

      /// Writes the xml file to the ostream 'os'
      void writeTo(std::ostream & os) override;

      // Docu in base class
      void onEndElement(const char16_t* qname) override;

      // Docu in base class
      void onStartElement(const char16_t* qname, const XMLAttributes& attributes) override;

      // Docu in base class
      void onCharacters(const char16_t* chars, Size /*length*/) override;

protected:
      std::map<std::string, std::pair<std::string, std::string> > & ptm_informations_;
      std::string name_, tag_, composition_;
      bool open_tag_;
    };

  }   // namespace Internal

} // namespace OpenMS

