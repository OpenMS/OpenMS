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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_RESIDUEDB_H
#define OPENMS_CHEMISTRY_RESIDUEDB_H

#include <OpenMS/DATASTRUCTURES/Map.h>
#include <boost/unordered_map.hpp>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <set>

namespace OpenMS
{
  // forward declarations
  class ResidueModification;
  class Residue;

  /** @ingroup Chemistry

      @brief residue data base which holds residues

      The residues stored in this DB are defined in a
      XML file under data/CHEMISTRY/residues.xml

      By default no modified residues are stored in an instance. However, if one
      queries the instance with getModifiedResidue, a new modified residue is
      added.
  */
  class OPENMS_DLLAPI ResidueDB
  {
public:

    /** @name Typedefs
    */
    //@{
    typedef std::set<Residue*>::iterator ResidueIterator;
    typedef std::set<const Residue*>::const_iterator ResidueConstIterator;
    //@}

    /// this member function serves as a replacement of the constructor
    inline static ResidueDB* getInstance()
    {
      static ResidueDB* db_ = nullptr;
      if (db_ == nullptr)
      {
        db_ = new ResidueDB;
      }
      return db_;
    }

    /** @name Constructors and Destructors
    */
    //@{
    /// destructor
    virtual ~ResidueDB();
    //@}

    /** @name Accessors
    */
    //@{
    /// returns the number of residues stored
    Size getNumberOfResidues() const;

    /// returns the number of modified residues stored
    Size getNumberOfModifiedResidues() const;

    /// returns a pointer to the residue with name, 3 letter code or 1 letter code name
    const Residue* getResidue(const String& name) const;

    /// returns a pointer to the residue with 1 letter code name
    const Residue* getResidue(const unsigned char& one_letter_code) const;

    /**
       @brief Returns a pointer to a modified residue given a modification name

       The "base" residue is looked up in ModificationsDB using the modification name.
       The modified residue is added to the database if it doesn't exist yet.
    */
    const Residue* getModifiedResidue(const String& name);

    /**
       @brief Returns a pointer to a modified residue given a residue and a modification name

       The modified residue is added to the database if it doesn't exist yet.

       @throw Exception::IllegalArgument if the residue was not found
       @throw Exception::InvalidValue if no matching modification was found (via ModificationsDB::getModification)
    */
    const Residue* getModifiedResidue(const Residue* residue, const String& name);

    /**
       @brief returns a set of all residues stored in this residue db

       The possible residues are defined in share/OpenMS/CHEMISTRY/Residues.xml.
       At the moment the following sets are available:
       All - all residues stored in the file
       Natural20 - default 20 naturally occurring residues
       Natural19WithoutI - default natural amino acids, excluding isoleucine (isobaric to leucine)
       Natural19WithoutL - default natural amino acids, excluding leucine (isobaric to isoleucine)
       Natural19J - default natural amino acids,  (isobaric leucine/isoleucine are marked by 'J')
       AmbiguousWithoutX - all amino acids, including ambiguous ones: B (asparagine or aspartate), Z (glutamine or glutamate), J (isoleucine or leucine)
       Ambiguous - all amino acids including all ambiguous ones (X can be every other amino acid)
       AllNatural - naturally occurring residues, including selenocysteine (U)

       @throw Exception::ElementNotFound if the specified residue set is not defined
    */
    const std::set<const Residue*> getResidues(const String& residue_set = "All") const;

    /// returns all residue sets that are registered which this instance
    const std::set<String>& getResidueSets() const;

    /// sets the residues from given file
    void setResidues(const String& filename);

    /// adds a residue, i.e. a unknown residue, where only the weight is known
    void addResidue(const Residue& residue);
    //@}

    /** @name Predicates
    */
    //@{
    /// returns true if the db contains a residue with the given name
    bool hasResidue(const String& name) const;

    /// returns true if the db contains the residue of the given pointer
    bool hasResidue(const Residue* residue) const;
    //@}

    /** @name Iterators
    */
    //@{
    inline ResidueIterator beginResidue() { return residues_.begin(); }

    inline ResidueIterator endResidue() { return residues_.end(); }

    inline ResidueConstIterator beginResidue() const { return const_residues_.begin(); }

    inline ResidueConstIterator endResidue() const { return const_residues_.end(); }
    //@}

protected:

    /** @name Private Constructors
    */
    //@{
    /// default constructor
    ResidueDB();

    ///copy constructor
    ResidueDB(const ResidueDB& residue_db);
    //@}

    /** @name Assignment
*/
    //@{
    /// assignment operator
    ResidueDB& operator=(const ResidueDB& aa);
    //@}

    /**
       @brief reads residues from the given file

       @throw Exception::ParseError if the file cannot be parsed
    */
    void readResiduesFromFile_(const String& filename);

    /// parses a residue, given the key/value pairs from i.e. an XML file
    Residue* parseResidue_(Map<String, String>& values);

    /// deletes all sub-instances of the stored data like modifications and residues
    void clear_();

    /// clears the residues
    void clearResidues_();

    /// builds an index of residue names for fast access, synonyms are also considered
    void buildResidueNames_();

    void addResidue_(Residue* residue);

    boost::unordered_map<String, Residue*> residue_names_;

    // fast lookup table for residues
    Residue* residue_by_one_letter_code_[256];

    Map<String, Map<String, Residue*> > residue_mod_names_;

    std::set<Residue*> residues_;

    std::set<const Residue*> const_residues_;

    std::set<Residue*> modified_residues_;

    std::set<const Residue*> const_modified_residues_;

    Map<String, std::set<const Residue*> > residues_by_set_;

    std::set<String> residue_sets_;
  };
}
#endif
