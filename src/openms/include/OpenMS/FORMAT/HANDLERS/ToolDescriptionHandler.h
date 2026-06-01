// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/DATASTRUCTURES/ToolDescription.h>
#include <OpenMS/FORMAT/HANDLERS/ParamXMLHandler.h>

namespace OpenMS
{
  class ProgressLogger;

  namespace Internal
  {

    /**
        @brief XML handler for ToolDescriptionFile

        @note Do not use this class. It is only needed in ToolDescriptionFile.
    */
    class OPENMS_DLLAPI ToolDescriptionHandler :
      public ParamXMLHandler
    {
public:
      /**@name Constructors and destructor */
      //@{

      /// Constructor
      ToolDescriptionHandler(const String & filename, const String & version);

      /// Destructor
      ~ToolDescriptionHandler() override;
      //@}


      // Docu in base class
      void endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname) override;

      // Docu in base class
      void startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes) override;

      // Docu in base class
      void characters(const XMLCh * const chars, const XMLSize_t length) override;

      // NOT IMPLEMENTED
      void writeTo(std::ostream & os) override;

      // Retrieve parsed tool description
      const std::vector<ToolDescription> & getToolDescriptions() const;

      // Set tool description for writing
      void setToolDescriptions(const std::vector<ToolDescription> & td);

protected:

      Param p_;

      Internal::ToolExternalDetails tde_;
      Internal::ToolDescription td_;
      std::vector<Internal::ToolDescription> td_vec_;

      String tag_;

      bool in_ini_section_;

private:

      ToolDescriptionHandler();
      ToolDescriptionHandler(const ToolDescriptionHandler & rhs);
      ToolDescriptionHandler & operator=(const ToolDescriptionHandler & rhs);

    };
  }   // namespace Internal
} // namespace OpenMS

