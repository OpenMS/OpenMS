// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer:  Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_TOOLDESCRIPTIONHANDLER_H
#define OPENMS_FORMAT_HANDLERS_TOOLDESCRIPTIONHANDLER_H

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

#endif // OPENMS_FORMAT_HANDLERS_TOOLDESCRIPTIONHANDLER_H
