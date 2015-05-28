// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Xiao Liang  $
// $Authors: Xiao Liang $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/EnzymesDB.h>
#include <OpenMS/CHEMISTRY/Enzyme.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>


#include <vector>
#include <fstream>

using namespace std;

namespace OpenMS
{
  EnzymesDB::EnzymesDB()
  {
    readEnzymesFromFile_("CHEMISTRY/Enzymes.xml");
    buildEnzymeNames_();
  }

  EnzymesDB::~EnzymesDB()
  {
    clear_();
  }

  const Enzyme* EnzymesDB::getEnzyme(const String& name) const
  {
    if (enzyme_names_.find(name) == enzyme_names_.end())
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Enzyme name cannot be found. '");
    }
    return enzyme_names_.at(name);
  }

  const Enzyme* EnzymesDB::getEnzymeByRegEx(const String & cleavage_regex) const
  {
    if (!enzyme_regex_.has(cleavage_regex))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Enzyme with regex "
                                                                                       + cleavage_regex + " was not registered in Enzyme DB, register first!").c_str());
    }
    return enzyme_regex_[cleavage_regex];
  }

  void EnzymesDB::setEnzymes(const String& file_name)
  {
    clearEnzymes_();
    readEnzymesFromFile_(file_name);
    buildEnzymeNames_();
  }

  void EnzymesDB::addEnzyme(const Enzyme& enzyme)
  {
    Enzyme* r = new Enzyme(enzyme);
    addEnzyme_(r);
  }

  void EnzymesDB::addEnzyme_(Enzyme* r)
  {
    enzymes_.insert(r);
    const_enzymes_.insert(r);
    vector<String> names;
    if (r->getName() != "")
    {
      names.push_back(r->getName());
    }
    set<String> synonyms = r->getSynonyms();
    for (set<String>::iterator it = synonyms.begin(); it != synonyms.end(); ++it)
    {
      names.push_back(*it);
    }
    buildEnzymeNames_();
    return;
  }

  bool EnzymesDB::hasRegEx(const String & cleavage_regex) const
  {
    if (enzyme_regex_.has(cleavage_regex))
    {
      return true;
    }
    return false;
  }

  bool EnzymesDB::hasEnzyme(const String& enzy_name) const
  {
    return (enzyme_names_.find(enzy_name) != enzyme_names_.end());
  }
  
  bool EnzymesDB::hasEnzyme(const Enzyme* enzyme) const
  {
    return (const_enzymes_.find(enzyme) != const_enzymes_.end() );
  }

  void EnzymesDB::readEnzymesFromFile_(const String& file_name)
  {
    String file = File::find(file_name);

    Param param;
    ParamXMLFile paramFile;
    paramFile.load(file, param);

    if (!param.begin().getName().hasPrefix("Enzymes"))
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", "");
    }

    try
    {
      vector<String> split;
      param.begin().getName().split(':', split);
      String prefix = split[0] + split[1];
      Enzyme* enzy_ptr = 0;
    
      Map<String, String> values;

      for (Param::ParamIterator it = param.begin(); it != param.end(); ++it)
      {
        it.getName().split(':', split);
        if (prefix != split[0] + split[1])
        {
          // add enzyme
          enzy_ptr = parseEnzyme_(values);
          values.clear();
          enzymes_.insert(enzy_ptr);
          const_enzymes_.insert(enzy_ptr);
          prefix = split[0] + split[1];
        }
        String value = it->value;
        String key = it.getName();
        values[key] = value;
      }

      // add last enzyme
      enzy_ptr = parseEnzyme_(values);
      enzymes_.insert(enzy_ptr);
      const_enzymes_.insert(enzy_ptr);
    }
    catch (Exception::BaseException& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, e.what(), "");
    }
  }

  void EnzymesDB::clear_()
  {
    clearEnzymes_();
  }

  void EnzymesDB::clearEnzymes_()
  {
    set<Enzyme*>::iterator it;
    for (it = enzymes_.begin(); it != enzymes_.end(); ++it)
    {
      delete *it;
    }
    enzymes_.clear();
    enzyme_names_.clear();
    enzyme_regex_.clear();
    const_enzymes_.clear();
  }

  Enzyme* EnzymesDB::parseEnzyme_(Map<String, String>& values)
  {
    Enzyme* enzy_ptr = new Enzyme("unknown_enzyme", "");

    for (Map<String, String>::iterator it = values.begin(); it != values.end(); ++it)
    {
      String key(it->first);
      String value(it->second);

      if (key.hasSuffix(":Name"))
      {
        enzy_ptr->setName(value);
        continue;
      }
      if (key.hasSuffix(":RegEx"))
      {
        enzy_ptr->setRegEx(value);
        continue;
      }
      if (key.hasSuffix(":RegExDescription"))
      {
        enzy_ptr->setRegExDescription(value);
        continue;
      }
      if (key.hasSuffix(":NTermGain"))
     {
        enzy_ptr->setNTermGain(EmpiricalFormula(value));
        continue;
      }
      if (key.hasSuffix(":CTermGain"))
      {
        enzy_ptr->setCTermGain(EmpiricalFormula(value));
        continue;
      }
      if (key.hasSubstring("PSIid"))
      {
        // no PSIid defined?
        if (!key.hasSuffix(":"))
        {
          enzy_ptr->setPSIid(value);
        }
        continue;
      }
      if (key.hasSubstring("XTANDEMid"))
      {
        if (!key.hasSuffix(":"))
        {
          enzy_ptr->setXTANDEMid(value);
        }
        continue;
      }
      if (key.hasSubstring("OMSSAid"))
      {
        if (!key.hasSuffix(":"))
        {
          enzy_ptr->setOMSSAid(value.toInt());
        }
        continue;
      }
      if (key.hasSubstring("Synonyms"))
      {
        // no synonyms defined?
        if (!key.hasSuffix(":"))
        {
          enzy_ptr->addSynonym(value);
        }
        continue;
      }
      cerr << "unknown key: " << key << ", with value: " << value << endl;
    }
    return enzy_ptr;
  }

  void EnzymesDB::buildEnzymeNames_()
  {
    set<Enzyme*>::iterator it;
    for (it = enzymes_.begin(); it != enzymes_.end(); ++it)
    {
      enzyme_names_[(*it)->getName()] = *it;
      if ((*it)->getRegEx() != "")
      {
        enzyme_regex_[(*it)->getRegEx()] = *it;
      }
    }
  }

  void EnzymesDB::getAllNames(vector<String> & all_names) const
  {
    all_names.clear();
    for (set<Enzyme *>::const_iterator it = enzymes_.begin(); it != enzymes_.end(); ++it)
    {
      all_names.push_back((*it)->getName());
    }
  }

  void EnzymesDB::getAllXTandemNames(vector<String> & all_names) const
  {
    all_names.clear();
    for (set<Enzyme *>::const_iterator it = enzymes_.begin(); it != enzymes_.end(); ++it)
    {
      if ((*it)->getXTANDEMid() != "")
      {
        all_names.push_back((*it)->getName());
      }
    }
  }

  void EnzymesDB::getAllOMSSANames(vector<String> & all_names) const
  {
    all_names.clear();
    all_names.push_back("Trypsin");
    for (set<Enzyme *>::const_iterator it = enzymes_.begin(); it != enzymes_.end(); ++it)
    {
      if ((*it)->getOMSSAid() != 0)
      {
        all_names.push_back((*it)->getName());
      }
    }
  }
} // namespace OpenMS
