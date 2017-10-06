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
// $Maintainer: Xiao Liang $
// $Authors: Xiao Liang, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_DIGESTIONENZYMEDB_H
#define OPENMS_CHEMISTRY_DIGESTIONENZYMEDB_H

#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <boost/unordered_map.hpp>
#include <set>

namespace OpenMS
{
  /**
    @ingroup Chemistry

    @brief Digestion enzyme database (base class)

    Template parameters:
    @p DigestionEnzymeType should be a subclass of DigestionEnzyme.
    @p InstanceType should be a subclass of DigestionEnzymeDB ("Curiously Recurring Template Pattern", see https://stackoverflow.com/a/34519373).
  */
  template<typename DigestionEnzymeType, typename InstanceType> class OPENMS_DLLAPI DigestionEnzymeDB
  {
  public:

    /** @name Typedefs
    */
    //@{
    typedef typename std::set<const DigestionEnzymeType*>::const_iterator ConstEnzymeIterator;
    typedef typename std::set<const DigestionEnzymeType*>::iterator EnzymeIterator;
    //@}

    /// this member function serves as a replacement of the constructor
    static InstanceType* getInstance()
    {
      static InstanceType* db_ = 0;
      if (db_ == 0)
      {
        db_ = new InstanceType;
      }
      return db_;
    }

    /** @name Constructors and Destructors
    */
    //@{
    /// destructor
    virtual ~DigestionEnzymeDB()
    {
      for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
      {
        delete *it;
      }
    }
    //@}

    /** @name Accessors
    */
    //@{
    /// returns a pointer to the enzyme with name (supports synonym names)
    /// @throw Exception::ElementNotFound if enzyme is unknown
    /// @note enzymes are registered in regular and in toLowercase() style, if unsure use toLowercase
    const DigestionEnzymeType* getEnzyme(const String& name) const
    {
      if (enzyme_names_.find(name) == enzyme_names_.end())
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
      }
      return enzyme_names_.at(name);
    }

    /// returns a pointer to the enzyme with cleavage regex
    /// @throw Exception::IllegalArgument if enzyme regex  is unregistered.
    const DigestionEnzymeType* getEnzymeByRegEx(const String& cleavage_regex) const
    {
      if (!enzyme_regex_.has(cleavage_regex))
      {
        // @TODO: why does this use a different exception than "getEnzyme"?
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         String("Enzyme with regex " + cleavage_regex + " was not registered in Enzyme DB, register first!").c_str());
      }
      return enzyme_regex_[cleavage_regex];
    }

    /// returns all the enzyme names (does NOT include synonym names)
    void getAllNames(std::vector<String>& all_names) const
    {
      all_names.clear();
      for (ConstEnzymeIterator it = const_enzymes_.begin(); it != const_enzymes_.end(); ++it)
      {
        all_names.push_back((*it)->getName());
      }
    }
    //@}

    /** @name Predicates
    */
    //@{
    /// returns true if the db contains a enzyme with the given name (supports synonym names)
    bool hasEnzyme(const String& name) const
    {
      return (enzyme_names_.find(name) != enzyme_names_.end());
    }

    /// returns true if the db contains a enzyme with the given regex
    bool hasRegEx(const String& cleavage_regex) const
    {
      return enzyme_regex_.has(cleavage_regex);
    }

    /// returns true if the db contains the enzyme of the given pointer
    bool hasEnzyme(const DigestionEnzymeType* enzyme) const
    {
      return (const_enzymes_.find(enzyme) != const_enzymes_.end() );
    }
    //@}

    /** @name Iterators
    */
    //@{
    inline ConstEnzymeIterator beginEnzyme() const { return const_enzymes_.begin(); }  // we only allow constant iterators -- this DB is not meant to be modifiable
    inline ConstEnzymeIterator endEnzyme() const { return const_enzymes_.end(); }

    //@}
  protected:
    DigestionEnzymeDB(const String& db_file = "")
    {
      if (!db_file.empty())
      {
        readEnzymesFromFile_(db_file);
      }
    }

    ///copy constructor
    DigestionEnzymeDB(const DigestionEnzymeDB& enzymes_db) = delete;
    //@}

    /** @name Assignment
    */
    //@{
    /// assignment operator
    DigestionEnzymeDB& operator=(const DigestionEnzymeDB& enzymes_db) = delete;
    //@}

    /// reads enzymes from the given file
    void readEnzymesFromFile_(const String& filename)
    {
      String file = File::find(filename);

      Param param;
      ParamXMLFile().load(file, param);
      if (param.empty()) return;

      std::vector<String> split;
      param.begin().getName().split(':', split);
      if (split[0] != "Enzymes")
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, split[0], "name 'Enzymes' expected");
      }

      try
      {
        Map<String, String> values;
        String previous_enzyme = split[1];
        // this iterates over all the "ITEM" elements in the XML file:
        for (Param::ParamIterator it = param.begin(); it != param.end(); ++it)
        {
          it.getName().split(':', split);
          if (split[0] != "Enzymes") break; // unexpected content in the XML file
          if (split[1] != previous_enzyme)
          {
            // add enzyme and reset:
            addEnzyme_(parseEnzyme_(values));
            previous_enzyme = split[1];
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

    /// parses an enzyme, given the key/value pairs from an XML file
    const DigestionEnzymeType* parseEnzyme_(Map<String, String>& values) const
    {
      DigestionEnzymeType* enzy_ptr = new DigestionEnzymeType();

      for (Map<String, String>::iterator it = values.begin(); it != values.end(); ++it)
      {
        const String& key = it->first;
        const String& value = it->second;
        if (!enzy_ptr->setValueFromFile(key, value))
        {
          LOG_ERROR << "Error while parsing enzymes file: unknown key '" << key << "' with value '" << value << "'" << std::endl;
        }
      }
      return enzy_ptr;
    }

    // add to internal data; also update indices for search by name and regex
    void addEnzyme_(const DigestionEnzymeType* enzyme)
    {
      // add to internal storage
      const_enzymes_.insert(enzyme);
      // add to internal indices (by name and its synonyms)
      String name = enzyme->getName();
      enzyme_names_[name] = enzyme;
      enzyme_names_[name.toLower()] = enzyme;
      for (std::set<String>::const_iterator it = enzyme->getSynonyms().begin(); it != enzyme->getSynonyms().end(); ++it)
      {
        enzyme_names_[*it] = enzyme;
      }
      // ... and by regex
      if (enzyme->getRegEx() != "")
      {
        enzyme_regex_[enzyme->getRegEx()] = enzyme;
      }
      return;
    }

    boost::unordered_map<String, const DigestionEnzymeType*> enzyme_names_;  // index by names

    Map<String, const DigestionEnzymeType*> enzyme_regex_; // index by regex

    std::set<const DigestionEnzymeType*> const_enzymes_; // set of enzymes

  };
}

#endif
