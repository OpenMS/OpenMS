// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/PTMXMLHandler.h>

using namespace std;
using namespace xercesc;

namespace OpenMS
{
  namespace Internal
  {
    PTMXMLHandler::PTMXMLHandler(map<String, pair<String, String> > & ptm_informations, const String & filename) :
      XMLHandler(filename, ""),
      ptm_informations_(ptm_informations)
    {
    }

    PTMXMLHandler::~PTMXMLHandler()
    {
    }

    void PTMXMLHandler::writeTo(std::ostream & os)
    {
      os << "<PTMs>" << "\n";
      for (map<String, pair<String, String> >::const_iterator ptm_i = ptm_informations_.begin(); ptm_i != ptm_informations_.end(); ++ptm_i)
      {
        os << "\t<PTM>" << "\n";
        os << "\t\t<name>" << ptm_i->first << "</name>" << "\n";             // see header
        os << "\t\t<composition>" << ptm_i->second.first << "</composition>" << "\n";
        os << "\t\t<possible_amino_acids>" << ptm_i->second.second << "</possible_amino_acids>" << "\n";
        os << "\t</PTM>" << "\n";
      }
      os << "</PTMs>" << "\n";
    }

    void PTMXMLHandler::startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const Attributes & /*attributes*/)
    {
      tag_ = String(sm_.convert(qname)).trim();
      open_tag_ = true;
    }

    void PTMXMLHandler::endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const /*qname*/)
    {
//          tag_ = String(sm_.convert(qname)).trim();
      tag_ = "";
      open_tag_ = false;
    }

    void PTMXMLHandler::characters(const XMLCh * const chars, const XMLSize_t /*length*/)
    {
      if (open_tag_)
      {
        if (tag_ == "name")
        {
          name_ = String(sm_.convert(chars)).trim();
        }
        else if (tag_ == "composition")
        {
          composition_ = String(sm_.convert(chars)).trim();
        }
        else if (tag_ == "possible_amino_acids")
        {
          ptm_informations_[name_] = make_pair(composition_, String(sm_.convert(chars)).trim());
        }
      }
    }

  }   // namespace Internal
} // namespace OpenMS
