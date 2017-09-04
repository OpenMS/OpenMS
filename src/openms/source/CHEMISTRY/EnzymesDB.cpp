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
// $Maintainer: Xiao Liang  $
// $Authors: Xiao Liang, Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/EnzymesDB.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeRNA.h>
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
    readEnzymesFromFile_<DigestionEnzymeRNA>("CHEMISTRY/Enzymes_RNA.xml");
  }

  EnzymesDB::~EnzymesDB()
  {
    clear();
  }

  const DigestionEnzyme* EnzymesDB::getEnzyme(const String& name) const
  {
    if (enzyme_names_.find(name) == enzyme_names_.end())
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Enzyme name cannot be found. '");
    }
    return enzyme_names_.at(name);
  }

  const DigestionEnzyme* EnzymesDB::getEnzymeByRegEx(const String& cleavage_regex) const
  {
    if (!enzyme_regex_.has(cleavage_regex))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Enzyme with regex "
                                                                                       + cleavage_regex + " was not registered in Enzyme DB, register first!").c_str());
    }
    return enzyme_regex_[cleavage_regex];
  }

  template <typename DigestionEnzymeType>
  void EnzymesDB::setEnzymes(const String& file_name)
  {
    clear();
    readEnzymesFromFile_<DigestionEnzymeType>(file_name);
  }

  void EnzymesDB::addEnzyme_(const DigestionEnzyme* r)
  {
    // add to internal storage
    const_enzymes_.insert(r);
    // add to internal indices (by name and its synonyms)
    String name = r->getName();
    enzyme_names_[name] = r;
    enzyme_names_[name.toLower()] = r;
    for (set<String>::const_iterator it = r->getSynonyms().begin(); it != r->getSynonyms().end(); ++it)
    {
      enzyme_names_[*it] = r;
    }
    // ... and by regex
    if (r->getRegEx() != "")
    {
      enzyme_regex_[r->getRegEx()] = r;
    }
    return;
  }

  bool EnzymesDB::hasRegEx(const String& cleavage_regex) const
  {
    return enzyme_regex_.has(cleavage_regex);
  }

  bool EnzymesDB::hasEnzyme(const String& enzy_name) const
  {
    return (enzyme_names_.find(enzy_name) != enzyme_names_.end());
  }

  bool EnzymesDB::hasEnzyme(const DigestionEnzyme* enzyme) const
  {
    return (const_enzymes_.find(enzyme) != const_enzymes_.end() );
  }

  template <typename DigestionEnzymeType>
  void EnzymesDB::readEnzymesFromFile_(const String& file_name)
  {
    String file = File::find(file_name);

    Param param;
    ParamXMLFile paramFile;
    paramFile.load(file, param);

    if (!param.begin().getName().hasPrefix("Enzymes"))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, param.begin().getName(), "prefix 'Enzymes' expected");
    }

    try
    {
      vector<String> split;
      param.begin().getName().split(':', split);
      String prefix = split[0] + split[1];

      Map<String, String> values;

      for (Param::ParamIterator it = param.begin(); it != param.end(); ++it)
      {
        it.getName().split(':', split);
        if (prefix != split[0] + split[1])
        {
          // add enzyme
          addEnzyme_(parseEnzyme_<DigestionEnzymeType>(values));
          prefix = split[0] + split[1];
          values.clear();
        }
        values[it.getName()] = it->value;
      }

      // add last enzyme
      addEnzyme_(parseEnzyme_(values));
    }
    catch (Exception::BaseException& e)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, e.what(), "");
    }
  }

  void EnzymesDB::clear()
  {
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      delete *it;
    }
    enzyme_names_.clear();
    enzyme_regex_.clear();
    const_enzymes_.clear();
  }

  template <typename DigestionEnzymeType>
  const DigestionEnzyme* EnzymesDB::parseEnzyme_(Map<String, String>& values) const
  {
    DigestionEnzyme* enzy_ptr = new DigestionEnzymeType();

    for (Map<String, String>::iterator it = values.begin(); it != values.end(); ++it)
    {
      const String& key = it->first;
      const String& value = it->second;
      if (!enzy_ptr->setValueFromFile_(key, value))
      {
        LOG_ERROR << "Error while parsing enzymes file: unknown key '" << key << "' with value '" << value << "'" << endl;
      }
    }
    return enzy_ptr;
  }

  void EnzymesDB::getAllNames(vector<String>& all_names) const
  {
    all_names.clear();
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      all_names.push_back((*it)->getName());
    }
  }

  void EnzymesDB::getAllXTandemNames(vector<String>& all_names) const
  {
    all_names.clear();
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      if ((*it)->getSubstrate() != DigestionEnzyme::PROTEIN) continue;

      const DigestionEnzymeProtein* p = dynamic_cast<const DigestionEnzymeProtein*>(*it);
      if (p->getXTandemID() != "")
      {
        all_names.push_back(p->getName());
      }
    }
  }

  void EnzymesDB::getAllCometNames(vector<String>& all_names) const
  {
    all_names.clear();
    all_names.push_back("unspecific cleavage");
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      if ((*it)->getSubstrate() != DigestionEnzyme::PROTEIN) continue;

      const DigestionEnzymeProtein* p = dynamic_cast<const DigestionEnzymeProtein*>(*it);
      if (p->getCometID() != 0)
      {
        all_names.push_back(p->getName());
      }
    }
  }

  void EnzymesDB::getAllOMSSANames(vector<String>& all_names) const
  {
    all_names.clear();
    all_names.push_back("Trypsin");
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      if ((*it)->getSubstrate() != DigestionEnzyme::PROTEIN) continue;

      const DigestionEnzymeProtein* p = dynamic_cast<const DigestionEnzymeProtein*>(*it);
      if (p->getOMSSAID() != 0)
      {
        all_names.push_back(p->getName());
      }
    }
  }

  void EnzymesDB::getAllMSGFNames(vector<String>& all_names) const
  {
    all_names.clear();
    for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
    {
      if ((*it)->getSubstrate() != DigestionEnzyme::PROTEIN) continue;

      const DigestionEnzymeProtein* p = dynamic_cast<const DigestionEnzymeProtein*>(*it);
      if (p->getMSGFID() != -1) // msgf+ starts enzyme numbering at 0
      {
        all_names.push_back(p->getName());
      }
    }
  }
} // namespace OpenMS
