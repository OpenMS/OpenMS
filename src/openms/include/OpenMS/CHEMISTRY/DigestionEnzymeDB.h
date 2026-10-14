// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Xiao Liang $
// $Authors: Xiao Liang, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeDataProvider.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <algorithm>
#include <cctype>
#include <set>
#include <map>
#include <memory>

namespace OpenMS
{
  /**
    @ingroup Chemistry

    @brief Digestion enzyme database (base class)

    Template parameters:
    @p DigestionEnzymeType should be a subclass of DigestionEnzyme.
    @p InstanceType should be a subclass of DigestionEnzymeDB ("Curiously Recurring Template Pattern", see https://stackoverflow.com/a/34519373).
  */
  template<typename DigestionEnzymeType, typename InstanceType> class DigestionEnzymeDB
  {
  public:

    /** @name Typedefs
    */
    //@{
    typedef typename std::set<const DigestionEnzymeType*>::const_iterator ConstEnzymeIterator;
    typedef typename std::set<const DigestionEnzymeType*>::iterator EnzymeIterator;
    //@}

    /// this member function serves as a replacement of the constructor
    /// The database is immutable after construction, so a const pointer is returned.
    static const InstanceType* getInstance()
    {
      static InstanceType* db_ = nullptr;
      if (db_ == nullptr)
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
    /// @note Names and synonyms are matched as written first, then ignoring case, whitespace, '-' and '_'
    ///       (e.g. "Glu-C", "glu_c" and "GluC" all find glutamyl endopeptidase). Other characters are significant,
    ///       as e.g. '/' and '+' change the specificity ("Trypsin/P" is not "Trypsin").
    const DigestionEnzymeType* getEnzyme(const std::string& name) const
    {
      const DigestionEnzymeType* enzyme = findEnzyme_(name);
      if (enzyme == nullptr)
      {
        throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, name);
      }
      return enzyme;
    }

    /// returns a pointer to the enzyme with cleavage regex
    /// @throw Exception::IllegalArgument if enzyme regex  is unregistered.
    const DigestionEnzymeType* getEnzymeByRegEx(const std::string& cleavage_regex) const
    {
      if (!hasRegEx(cleavage_regex))
      {
        // @TODO: why does this use a different exception than "getEnzyme"?
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         std::string("Enzyme with regex " + cleavage_regex + " was not registered in Enzyme DB, register first!").c_str());
      }
      return enzyme_regex_.at(cleavage_regex);
    }

    /// returns all the enzyme names (does NOT include synonym names)
    void getAllNames(std::vector<std::string>& all_names) const
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
    /// returns true if the db contains a enzyme with the given name (supports synonym names; matched as in getEnzyme())
    bool hasEnzyme(const std::string& name) const
    {
      return (findEnzyme_(name) != nullptr);
    }

    /// returns true if the db contains a enzyme with the given regex
    bool hasRegEx(const std::string& cleavage_regex) const
    {
      return (enzyme_regex_.contains(cleavage_regex));
    }

    /// returns true if the db contains the enzyme of the given pointer
    bool hasEnzyme(const DigestionEnzymeType* enzyme) const
    {
      return (const_enzymes_.contains(enzyme) );
    }
    //@}

    /** @name Iterators
    */
    //@{
    inline ConstEnzymeIterator beginEnzyme() const { return const_enzymes_.begin(); }  // we only allow constant iterators -- this DB is not meant to be modifiable
    inline ConstEnzymeIterator endEnzyme() const { return const_enzymes_.end(); }

    //@}
  protected:
    DigestionEnzymeDB() = default;

    ///copy constructor
    DigestionEnzymeDB(const DigestionEnzymeDB& enzymes_db) = delete;
    //@}

    /** @name Assignment
    */
    //@{
    /// assignment operator
    DigestionEnzymeDB& operator=(const DigestionEnzymeDB& enzymes_db) = delete;
    //@}

    /// add to internal data; also update indices for search by name and regex.
    /// If an enzyme with the same name already exists, it is replaced.
    void addEnzyme_(const DigestionEnzymeType* enzyme)
    {
      std::string name = enzyme->getName();

      // if an enzyme with the same name exists, remove the old one first
      // (compare names only: a synonym or lookup key of another enzyme must not remove that enzyme)
      auto existing = std::find_if(const_enzymes_.begin(), const_enzymes_.end(),
                                   [&name](const DigestionEnzymeType* e) { return e->getName() == name; });
      if (existing != const_enzymes_.end())
      {
        const DigestionEnzymeType* old = *existing;
        const_enzymes_.erase(existing);
        // remove old name/synonym entries, but only those that point to the enzyme being replaced
        auto points_to_old = [old](const auto& entry) { return entry.second == old; };
        std::erase_if(enzyme_names_, points_to_old);
        std::erase_if(enzyme_keys_, points_to_old);
        // remove old regex entry
        if (!old->getRegEx().empty())
        {
          enzyme_regex_.erase(old->getRegEx());
        }
        delete old;
      }

      // add to internal storage
      const_enzymes_.insert(enzyme);
      // add to internal indices (by name and its synonyms, as written ...
      enzyme_names_[name] = enzyme;
      for (const auto& syn : enzyme->getSynonyms())
      {
        enzyme_names_[syn] = enzyme;
      }
      // ... and as lookup keys, see getEnzyme(); a synonym never re-points an existing key)
      enzyme_keys_[lookupKey_(name)] = enzyme;
      for (const auto& syn : enzyme->getSynonyms())
      {
        enzyme_keys_.emplace(lookupKey_(syn), enzyme);
      }
      // ... and by regex
      if (enzyme->getRegEx() != "")
      {
        enzyme_regex_[enzyme->getRegEx()] = enzyme;
      }
      return;
    }

    /// Load enzymes from a list of data providers and add them to the database.
    /// Each provider's loadEnzymes() is called in order; enzymes with duplicate names
    /// from later providers replace earlier ones (allowing user overrides).
    void loadFromProviders_(std::vector<std::unique_ptr<DigestionEnzymeDataProvider<DigestionEnzymeType>>>& providers)
    {
      for (auto& provider : providers)
      {
        auto enzymes = provider->loadEnzymes();
        for (auto& enzyme : enzymes)
        {
          addEnzyme_(enzyme.release());
        }
      }
    }

    /// key for lenient name lookup (see getEnzyme()): lowercase, without whitespace, '-' and '_'
    static std::string lookupKey_(const std::string& name)
    {
      std::string key;
      key.reserve(name.size());
      for (unsigned char c : name)
      {
        if (std::isspace(c) || c == '-' || c == '_') continue;
        key += static_cast<char>(std::tolower(c));
      }
      return key;
    }

    /// returns the enzyme with the given name or synonym (see getEnzyme()), or nullptr if unknown
    const DigestionEnzymeType* findEnzyme_(const std::string& name) const
    {
      auto pos = enzyme_names_.find(name);
      if (pos != enzyme_names_.end()) return pos->second;
      pos = enzyme_keys_.find(lookupKey_(name));
      return (pos != enzyme_keys_.end()) ? pos->second : nullptr;
    }

    std::map<std::string, const DigestionEnzymeType*> enzyme_names_; ///< index by names and synonyms (as written)

    std::map<std::string, const DigestionEnzymeType*> enzyme_keys_; ///< index by lookupKey_() of names and synonyms

    std::map<std::string, const DigestionEnzymeType*> enzyme_regex_; ///< index by regex

    std::set<const DigestionEnzymeType*> const_enzymes_; ///< set of enzymes

  };
}

