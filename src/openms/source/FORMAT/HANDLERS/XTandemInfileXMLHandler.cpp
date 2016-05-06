// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/XTandemInfileXMLHandler.h>
#include <xercesc/sax2/Attributes.hpp>

using namespace std;
using namespace xercesc;

namespace OpenMS
{
  namespace Internal
  {

    XTandemInfileXMLHandler::XTandemInfileXMLHandler(const String & filename, vector<XTandemInfileNote> & notes) :
      XMLHandler(filename, ""),
      notes_(notes), // this is a reference!
      actual_note_(),
      tag_()
    {
    }

    XTandemInfileXMLHandler::~XTandemInfileXMLHandler()
    {
    }

    void XTandemInfileXMLHandler::startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const Attributes & attributes)
    {

      tag_.push_back(String(sm_.convert(qname)));

      if (tag_.back() == "note")
      {
        int type_idx = attributes.getIndex(sm_.convert("type"));
        int label_idx = attributes.getIndex(sm_.convert("label"));

        if (type_idx != -1)
        {
          actual_note_.note_type = String(sm_.convert(attributes.getValue(type_idx)));
        }
        if (label_idx != -1)
        {
          actual_note_.note_label = String(sm_.convert(attributes.getValue(label_idx)));
        }
      }

    }

    void XTandemInfileXMLHandler::endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname)
    {
      String tag_close = String(sm_.convert(qname)).trim();
      if (tag_.back() != tag_close)
      {
        fatalError(LOAD, "Invalid closing/opening tag sequence. Unexpected tag '</ " + tag_close + ">'!");
      }
      if (tag_.back() == "note")
      {
        notes_.push_back(actual_note_);
        // prepare for new note
        actual_note_ = XTandemInfileNote();
      }
      
      tag_.pop_back();
    }

    void XTandemInfileXMLHandler::characters(const XMLCh * const chars, const XMLSize_t /*length*/)
    {
      if (tag_.back() == "note")
      {
        actual_note_.note_value = ((String)sm_.convert(chars)).trim();
      }
    }

  }   // namespace Internal
} // namespace OpenMS
