// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// NOTE: This is a private (non-installed) libOpenMS header. It bridges the
// Xerces SAX2 interface to OpenMS' Xerces-free XMLHandler callbacks and must only
// be included from libOpenMS .cpp translation units (never from a public header).

#pragma once

#include <OpenMS/FORMAT/HANDLERS/StringManager.h>
#include <OpenMS/FORMAT/HANDLERS/XMLAttributes.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/sax/SAXParseException.hpp>

namespace OpenMS::Internal
{
  /**
    @brief Adapts an OpenMS @ref XMLHandler to the Xerces SAX2 @c DefaultHandler interface.

    Xerces owns the SAX callback interface; OpenMS handlers expose the Xerces-free
    @c onStartElement / @c onEndElement / @c onCharacters callbacks (and OpenMS-typed
    error handlers). This adapter is installed as the Xerces content+error handler and
    forwards each callback to the wrapped @ref XMLHandler. XMLCh is guaranteed to be
    char16_t-sized, so the reinterpret_casts at this boundary are well-defined views.

    The forwarding methods are intentionally not @c noexcept: subclass callbacks may
    throw (e.g. XMLHandler::EndParsingSoftly to abort parsing early), and that
    exception must unwind through Xerces' parse() exactly as before.
  */
  class SAX2HandlerAdapter : public xercesc::DefaultHandler
  {
  public:
    explicit SAX2HandlerAdapter(XMLHandler& handler) : handler_(handler) {}

    void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*localname*/,
                      const XMLCh* const qname, const xercesc::Attributes& attrs) override
    {
      handler_.onStartElement(reinterpret_cast<const char16_t*>(qname), XMLAttributes(attrs));
    }

    void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*localname*/,
                    const XMLCh* const qname) override
    {
      handler_.onEndElement(reinterpret_cast<const char16_t*>(qname));
    }

    void characters(const XMLCh* const chars, const XMLSize_t length) override
    {
      handler_.onCharacters(reinterpret_cast<const char16_t*>(chars), length);
    }

    void fatalError(const xercesc::SAXParseException& e) override
    {
      handler_.fatalError(XMLHandler::LOAD, message_(e), (UInt)e.getLineNumber(), (UInt)e.getColumnNumber());
    }

    void error(const xercesc::SAXParseException& e) override
    {
      handler_.error(XMLHandler::LOAD, message_(e), (UInt)e.getLineNumber(), (UInt)e.getColumnNumber());
    }

    void warning(const xercesc::SAXParseException& e) override
    {
      handler_.warning(XMLHandler::LOAD, message_(e), (UInt)e.getLineNumber(), (UInt)e.getColumnNumber());
    }

  private:
    static std::string message_(const xercesc::SAXParseException& e)
    {
      return StringManager::convert(reinterpret_cast<const char16_t*>(e.getMessage()));
    }

    XMLHandler& handler_;
  };

} // namespace OpenMS::Internal
