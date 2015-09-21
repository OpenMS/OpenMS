
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
// $Maintainer: Xiao Liang $
// $Authors: Xiao Liang, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_ENZYMESDB_H
#define OPENMS_CHEMISTRY_ENZYMESDB_H

#include <OpenMS/DATASTRUCTURES/Map.h>
#include <boost/unordered_map.hpp>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <set>

namespace OpenMS
{
  // forward declarations
  class Enzyme;

  /** @ingroup Chemistry

    @brief enzyme database which holds enzymes

    The enzymes stored in this DB are defined in a
    XML file under share/CHEMISTRY/Enzymes.xml

  */
  class OPENMS_DLLAPI EnzymesDB
  {
public:

    /** @name Typedefs
    */
    //@{
    typedef std::set<const Enzyme *>::const_iterator ConstEnzymeIterator;
    typedef std::set<const Enzyme *>::iterator EnzymeIterator;
    //@}

    /// this member function serves as a replacement of the constructor
    inline static EnzymesDB* getInstance()
    {
      static EnzymesDB * db_ = 0;
      if (db_ == 0)
      {
        db_ = new EnzymesDB;
      }
      return db_;
    }

    /** @name Constructors and Destructors
    */
    //@{
    /// destructor
    virtual ~EnzymesDB();
    //@}

    /** @name Accessors
    */
    //@{
    /// returns a pointer to the enzyme with name (supports synonym names)
    /// @throw Exception::ElementNotFound if enzyme is unknown
    const Enzyme* getEnzyme(const String & name) const;

    /// returns a pointer to the enzyme with cleavage regex
    /// @throw Exception::IllegalArgument if enzyme regex  is unregistered.
    const Enzyme* getEnzymeByRegEx(const String & cleavage_regex) const;

    /// load enzymes from given file
    void setEnzymes(const String& filename);

    /// adds a new enzyme, e.g. a new enzyme, where only the cleavage regex is known
    void addEnzyme(const Enzyme& enzyme);

    /// deletes all enzymes, resulting in an empty database
    void clear();

    /// returns all the enzyme names (does NOT include synonym names)
    void getAllNames(std::vector<String>& all_names) const;

    /// returns all the enzyme names available for XTandem
    void getAllXTandemNames(std::vector<String>& all_names) const;
    
    /// returns all the enzyme names available for OMSSA
    void getAllOMSSANames(std::vector<String>& all_names) const;


    //@}


    /** @name Predicates
    */
    //@{
    /// returns true if the db contains a enzyme with the given name (supports synonym names)
    bool hasEnzyme(const String& name) const;

    /// returns true if the db contains a enzyme with the given regex
    bool hasRegEx(const String& cleavage_regex) const;
    
    /// returns true if the db contains the enzyme of the given pointer
    bool hasEnzyme(const Enzyme* enzyme) const;
    //@}

    /** @name Iterators
    */
    //@{
    inline ConstEnzymeIterator beginEnzyme() const { return const_enzymes_.begin(); }  // we only allow constant iterators -- this DB is not meant to be modifiable
    inline ConstEnzymeIterator endEnzyme() const { return const_enzymes_.end(); }

    //@}
protected:
    EnzymesDB();
    
    ///copy constructor
    EnzymesDB(const EnzymesDB& enzyme_db);
    //@}

    /** @name Assignment
    */
    //@{
    /// assignment operator
    EnzymesDB & operator=(const EnzymesDB& enzymes_db);
    //@}

    /// reads enzymes from the given file
    void readEnzymesFromFile_(const String& filename);

    /// parses a enzyme, given the key/value pairs from i.e. an XML file
    const Enzyme* parseEnzyme_(Map<String, String>& values) const;

    // add to internal data; also update indices for search by name and regex
    void addEnzyme_(const Enzyme* enzyme);

    boost::unordered_map<String, const Enzyme*> enzyme_names_;  // index by names

    Map<String, const Enzyme*> enzyme_regex_; // index by regex

    std::set<const Enzyme*> const_enzymes_; // set of enzymes

  };
}
#endif
