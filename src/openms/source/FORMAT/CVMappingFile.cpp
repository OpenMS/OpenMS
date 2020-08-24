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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/DATASTRUCTURES/CVReference.h>
#include <OpenMS/DATASTRUCTURES/CVMappingTerm.h>
#include <OpenMS/SYSTEM/File.h>

using namespace xercesc;
using namespace std;

namespace OpenMS
{

  CVMappingFile::CVMappingFile() :
    XMLHandler("", 0),
    XMLFile()
  {

  }

  CVMappingFile::~CVMappingFile()
  {
  }

  void CVMappingFile::load(const String& filename, CVMappings& cv_mappings, bool strip_namespaces)
  {
    //File name for error messages in XMLHandler
    file_ = filename;

    strip_namespaces_ = strip_namespaces;

    parse_(filename, this);

    cv_mappings.setCVReferences(cv_references_);
    cv_mappings.setMappingRules(rules_);

    cv_references_.clear();
    rules_.clear();

    return;
  }

  void CVMappingFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
  {

    tag_ = String(sm_.convert(qname));

    if (tag_ == "CvReference")
    {
      // CvReference cvName="PSI-PI" cvIdentifier="PSI-PI"/>
      CVReference ref;
      ref.setName(attributeAsString_(attributes, "cvName"));
      ref.setIdentifier(attributeAsString_(attributes, "cvIdentifier"));
      cv_references_.push_back(ref);
      return;
    }

    if (tag_ == "CvMappingRule")
    {
      // id="R1" cvElementPath="/psi-pi:MzIdentML/psi-pi:AnalysisSoftwareList/psi-pi:AnalysisSoftware/pf:ContactRole/pf:role/pf:cvParam" requirementLevel="MUST"  scopePath="" cvTermsCombinationLogic="OR
      actual_rule_.setIdentifier(attributeAsString_(attributes, "id"));
      String element_path = attributeAsString_(attributes, "cvElementPath");
      if (strip_namespaces_)
      {
        vector<String> slash_split;
        element_path.split('/', slash_split);
        if (slash_split.empty())
        {
          slash_split.push_back(element_path);
        }
        element_path = "";
        for (vector<String>::const_iterator it = slash_split.begin(); it != slash_split.end(); ++it)
        {
          if (it->empty())
          {
            continue;
          }

          vector<String> split;
          it->split(':', split);
          if (split.empty())
          {
            element_path += "/" + *it;
          }
          else
          {
            if (split.size() == 2)
            {
              element_path += "/" + split[1];
            }
            else
            {
              fatalError(LOAD, String("Cannot parse namespaces of path: '") + element_path + "'");
            }
          }
        }
      }
      actual_rule_.setElementPath(element_path);
      CVMappingRule::RequirementLevel level = CVMappingRule::MUST;
      String lvl = attributeAsString_(attributes, "requirementLevel");
      if (lvl == "MAY")
      {
        level = CVMappingRule::MAY;
      }
      else
      {
        if (lvl == "SHOULD")
        {
          level = CVMappingRule::SHOULD;
        }
        else
        {
          if (lvl == "MUST")
          {
            level = CVMappingRule::MUST;
          }
          else
          {
            // throw Exception
          }
        }
      }

      actual_rule_.setRequirementLevel(level);

      actual_rule_.setScopePath(attributeAsString_(attributes, "scopePath"));
      CVMappingRule::CombinationsLogic logic = CVMappingRule::OR;
      String lgc = attributeAsString_(attributes, "cvTermsCombinationLogic");
      if (lgc == "OR")
      {
        logic = CVMappingRule::OR;
      }
      else
      {
        if (lgc == "AND")
        {
          logic = CVMappingRule::AND;
        }
        else
        {
          if (lgc == "XOR")
          {
            logic = CVMappingRule::XOR;
          }
          else
          {
            // throw Exception;
          }
        }
      }
      actual_rule_.setCombinationsLogic(logic);
      return;
    }

    if (tag_ == "CvTerm")
    {
      // termAccession="PI:00266" useTermName="false" useTerm="false" termName="role type" isRepeatable="true" allowChildren="true" cvIdentifierRef="PSI-PI"
      CVMappingTerm term;

      term.setAccession(attributeAsString_(attributes, "termAccession"));
      term.setUseTerm(DataValue(attributeAsString_(attributes, "useTerm")).toBool());

      String use_term_name;
      optionalAttributeAsString_(use_term_name, attributes, "useTermName");
      if (use_term_name != "")
      {
        term.setUseTermName(DataValue(use_term_name).toBool());
      }
      else
      {
        term.setUseTermName(false);
      }
      term.setTermName(attributeAsString_(attributes, "termName"));

      String is_repeatable;
      optionalAttributeAsString_(is_repeatable, attributes, "isRepeatable");
      if (is_repeatable != "")
      {
        term.setIsRepeatable(DataValue(is_repeatable).toBool());
      }
      else
      {
        term.setIsRepeatable(true);
      }
      term.setAllowChildren(DataValue(attributeAsString_(attributes, "allowChildren")).toBool());
      term.setCVIdentifierRef(attributeAsString_(attributes, "cvIdentifierRef"));

      actual_rule_.addCVTerm(term);
      return;
    }

    return;
  }

  void CVMappingFile::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
  {
    tag_ = String(sm_.convert(qname));

    if (tag_ == "CvMappingRule")
    {
      rules_.push_back(actual_rule_);
      actual_rule_ = CVMappingRule();
      return;
    }

    return;
  }

  void CVMappingFile::characters(const XMLCh* const /*chars*/, const XMLSize_t /*length*/)
  {
    // good XML format, nothing to do here
  }

} // namespace OpenMS
