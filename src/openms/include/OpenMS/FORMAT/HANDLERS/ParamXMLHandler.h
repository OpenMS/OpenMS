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
      ParamXMLHandler(Param& param, const std::string& filename, const std::string& version);
      /// Destructor
      ~ParamXMLHandler() override;

      // Docu in base class
      void onEndElement(const char16_t* qname) override;

      // Docu in base class
      void onStartElement(const char16_t* qname, const XMLAttributes& attributes) override;

      /**
        @brief Whether a @em PARAMETERS element was encountered while parsing.

        Callers use this to reject documents that are not parameter files at all: unknown elements
        are silently ignored, so without this check any well-formed XML parses "successfully" into
        an empty Param. Note that this deliberately does not test the *root* element -- CTD files
        wrap @em PARAMETERS in a @em tool root, and ToolDescriptionHandler embeds it in @em ttd.
      */
      bool hasSeenParametersElement() const
      {
        return parameters_element_seen_;
      }

protected:
      /// The current absolute path (concatenation of nodes_ with <i>:</i> in between)
      std::string path_;
      /// Reference to the Param object to fill
      Param& param_;
      /// Map of node descriptions (they are set at the end of parsing)
      std::map<std::string, std::string> descriptions_;

      /// Whether a 'PARAMETERS' element was seen (see hasSeenParametersElement())
      bool parameters_element_seen_{false};

      ///Temporary data for parsing of item lists
      struct
      {
        std::string name;
        std::string type;
        std::vector<std::string> stringlist;
        IntList intlist;
        DoubleList doublelist;
        std::vector<std::string> tags;
        std::string description;
        std::string restrictions;
        Int restrictions_index;
      } list_;

private:
      /// Not implemented
      ParamXMLHandler();
    };

  } // namespace Internal
} // namespace OpenMS

