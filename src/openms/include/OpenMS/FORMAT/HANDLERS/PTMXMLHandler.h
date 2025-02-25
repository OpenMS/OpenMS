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
      PTMXMLHandler(std::map<String, std::pair<String, String> > & ptm_informations, const String & filename);

      /// Destructor
      ~PTMXMLHandler() override;

      /// Writes the xml file to the ostream 'os'
      void writeTo(std::ostream & os) override;

      // Docu in base class
      void endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname) override;

      // Docu in base class
      void startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes) override;

      // Docu in base class
      void characters(const XMLCh * const chars, const XMLSize_t /*length*/) override;

protected:
      std::map<String, std::pair<String, String> > & ptm_informations_;
      String name_, tag_, composition_;
      bool open_tag_;
    };

  }   // namespace Internal

} // namespace OpenMS

