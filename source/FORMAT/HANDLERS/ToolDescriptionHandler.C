// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#include <OpenMS/FORMAT/HANDLERS/ToolDescriptionHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <set>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {

    ToolDescriptionHandler::ToolDescriptionHandler(const String & filename, const String & version) :
      ParamXMLHandler(p_, filename, version),
      p_(),
      tde_(),
      td_(),
      td_vec_(),
      tag_(),
      in_ini_section_(false)
    {
    }

    ToolDescriptionHandler::~ToolDescriptionHandler()
    {
    }

    void ToolDescriptionHandler::startElement(const XMLCh * const uri, const XMLCh * const local_name, const XMLCh * const qname, const xercesc::Attributes & attributes)
    {
      if (in_ini_section_)
      {
        ParamXMLHandler::startElement(uri, local_name, qname, attributes);
        return;
      }

      tag_ = sm_.convert(qname);
      open_tags_.push_back(tag_);
      //std::cout << "starting tag " << tag_ << "\n";

      if (tag_ == "tool")
      {
        String status = attributeAsString_(attributes, "status");
        if (status == "external")
          td_.is_internal = false;
        else if (status == "internal")
          td_.is_internal = true;
        else
          error(LOAD, "ToolDescriptionHandler::startElement: Element 'status' if tag 'tool' has unknown value " + status + "'.");
        return;
      }
      if (tag_ == "mapping")
      {
        Int id = attributeAsInt_(attributes, "id");
        String command = attributeAsString_(attributes, "cl");
        tde_.tr_table.mapping[id] = command;
        return;
      }
      if (tag_ == "file_post")
      {
        Internal::FileMapping fm;
        fm.location = attributeAsString_(attributes, "location");
        fm.target = attributeAsString_(attributes, "target");
        tde_.tr_table.post_moves.push_back(fm);
        return;
      }
      if (tag_ == "file_pre")
      {
        Internal::FileMapping fm;
        fm.location = attributeAsString_(attributes, "location");
        fm.target = attributeAsString_(attributes, "target");
        tde_.tr_table.pre_moves.push_back(fm);
        return;
      }

      if (tag_ == "ini_param")
      {
        in_ini_section_ = true;
        p_ = Param(); // reset Param
        return;
      }

      if (tag_ == "ttd" || tag_ == "category" || tag_ == "e_category" || tag_ == "type")
        return;

      if (td_.is_internal)
      {
        if (tag_ == "name")
        {
          return;
        }
      }
      else if (!td_.is_internal)
      {
        if (tag_ == "external" || tag_ == "cloptions" || tag_ == "path" || tag_ == "mappings" || tag_ == "mapping" || tag_ == "ini_param" ||
            tag_ == "text" || tag_ == "onstartup" || tag_ == "onfail" || tag_ == "onfinish" || tag_ == "workingdirectory")
          return;
      }

      error(LOAD, "ToolDescriptionHandler::startElement(): Unkown element found: '" + tag_ + "', ignoring.");
    }

    void ToolDescriptionHandler::characters(const XMLCh * const chars, const XMLSize_t length)
    {
      if (in_ini_section_)
      {
        ParamXMLHandler::characters(chars, length);
        return;
      }

      //std::cout << "characters '" << sm_.convert(chars) << "' in tag " << tag_ << "\n";

      if (tag_ == "ttd" || tag_ == "tool" || tag_ == "mappings" || tag_ == "external" || tag_ == "text")
        return;

      if (tag_ == "name")
        td_.name = sm_.convert(chars);
      else if (tag_ == "category")
        td_.category = sm_.convert(chars);
      else if (tag_ == "type")
        td_.types.push_back(sm_.convert(chars));
      else if (tag_ == "e_category")
        tde_.category = sm_.convert(chars);
      else if (tag_ == "cloptions")
        tde_.commandline = sm_.convert(chars);
      else if (tag_ == "path")
        tde_.path = sm_.convert(chars);
      else if (tag_ == "onstartup")
        tde_.text_startup = sm_.convert(chars);
      else if (tag_ == "onfail")
        tde_.text_fail = sm_.convert(chars);
      else if (tag_ == "onfinish")
        tde_.text_finish = sm_.convert(chars);
      else if (tag_ == "workingdirectory")
        tde_.working_directory = sm_.convert(chars);

      else
        error(LOAD, "ToolDescriptionHandler::characters: Unkown character section found: '" + tag_ + "', ignoring.");
    }

    void ToolDescriptionHandler::endElement(const XMLCh * const uri, const XMLCh * const local_name, const XMLCh * const qname)
    {
      String endtag_ = sm_.convert(qname);
      if (in_ini_section_ && endtag_ != "ini_param")
      {
        ParamXMLHandler::endElement(uri, local_name, qname);
        return;
      }

      open_tags_.pop_back();
      //std::cout << "ending tag " << endtag_ << "\n";
      if (open_tags_.size() > 0)
        tag_ = open_tags_.back();
      //std::cout << " --> current Tag: " << tag_ << "\n";

      if (endtag_ == "ini_param")
      {
        in_ini_section_ = false;
        tde_.param = p_;
        return;
      }
      else if (endtag_ == "external")
      {
        td_.external_details.push_back(tde_);
        tde_ = ToolExternalDetails();
        return;
      }
      else if (endtag_ == "tool")
      {
        td_vec_.push_back(td_);
        td_ = ToolDescription();
        return;
      }
      else
        return;  // TODO...

      error(LOAD, "ToolDescriptionHandler::endElement: Unkown element found: '" + endtag_ + "', ignoring.");

    }

    void ToolDescriptionHandler::writeTo(std::ostream & /*os*/)
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }

    void ToolDescriptionHandler::setToolDescriptions(const std::vector<ToolDescription> & tds)
    {
      td_vec_ = tds;
    }

    const std::vector<ToolDescription> & ToolDescriptionHandler::getToolDescriptions() const
    {
      return td_vec_;
    }

  }   //namespace Internal
} // namespace OpenMS
