// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Stephan Aiche $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

#include <vector>
#include <map>

namespace OpenMS
{
  namespace Internal
  {
    /**
        @brief XML Handler for Param files.

    */
    class OPENMS_DLLAPI ParamXMLHandler :
      public XMLHandler
    {
public:
      /// Default constructor
      ParamXMLHandler(Param& param, const String& filename, const String& version);
      /// Destructor
      ~ParamXMLHandler() override;

      // Docu in base class
      void endElement(const XMLCh* const uri, const XMLCh* const local_name, const XMLCh* const qname) override;

      // Docu in base class
      void startElement(const XMLCh* const uri, const XMLCh* const local_name, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

protected:
      /// The current absolute path (concatenation of nodes_ with <i>:</i> in between)
      String path_;
      /// Reference to the Param object to fill
      Param& param_;
      /// Map of node descriptions (they are set at the end of parsing)
      std::map<String, String> descriptions_;

      ///Temporary data for parsing of item lists
      struct
      {
        String name;
        String type;
        std::vector<std::string> stringlist;
        IntList intlist;
        DoubleList doublelist;
        std::vector<std::string> tags;
        String description;
        String restrictions;
        Int restrictions_index;
      } list_;

private:
      /// Not implemented
      ParamXMLHandler();
    };

  } // namespace Internal
} // namespace OpenMS

