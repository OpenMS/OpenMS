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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Stephan Aiche $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_PARAMXMLHANDLER_H
#define OPENMS_FORMAT_HANDLERS_PARAMXMLHANDLER_H

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
        StringList stringlist;
        IntList intlist;
        DoubleList doublelist;
        StringList tags;
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

#endif // OPENMS_FORMAT_HANDLERS_PARAMXMLHANDLER_H
