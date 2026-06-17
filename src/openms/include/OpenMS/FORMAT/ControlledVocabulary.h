// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Andreas Bertsch, Mathias Walzer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/ListUtils.h> // StringList
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <set>
#include <map>

namespace OpenMS
{
  /**
      @brief Representation of a controlled vocabulary.

      This representation only contains the information used for parsing and validation.
      All other lines are stored in the @em unparsed member of the CVTerm struct.

  @ingroup Format
  */
  class OPENMS_DLLAPI ControlledVocabulary
  {
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const ControlledVocabulary& cv);

public:
    /// ensure same hash on all platforms (for reproducibility)-
    struct FNV1aHasher
    {
      size_t operator()(const std::string& key) const noexcept
      {
        size_t hash = 14695981039346656037ull;
        for (auto c : key)
        {
          hash ^= static_cast<unsigned char>(c);
          hash *= 1099511628211ull;
        }
        return hash;
      }
    };

    /// Representation of a CV term
    struct OPENMS_DLLAPI CVTerm
    {
      /// define xsd types allowed in cv term to specify their value-type
      enum class XRefType
      {
        XSD_STRING = 0, // xsd:string A string
        XSD_INTEGER, // xsd:integer Any integer
        XSD_DECIMAL, // xsd:decimal Any real number
        XSD_NEGATIVE_INTEGER, // xsd:negativeInteger Any negative integer
        XSD_POSITIVE_INTEGER, // xsd:positiveInteger Any integer > 0
        XSD_NON_NEGATIVE_INTEGER, // xsd:nonNegativeInteger Any integer >= 0
        XSD_NON_POSITIVE_INTEGER, // xsd:nonPositiveInteger Any integer < 0
        XSD_BOOLEAN, // xsd:boolean True or false
        XSD_DATE, // xsd:date An XML-Schema date
        XSD_ANYURI, // xsd:anyURI uniform resource identifier
        NONE
      };

      static std::string getXRefTypeName(XRefType type);
      //static bool isSearchEngineSpecificScore();
      static bool isHigherBetterScore(ControlledVocabulary::CVTerm term); ///if it is a score type, lookup has_order

      std::string name; ///< Text name
      std::string id; ///< Identifier
      std::set<std::string> parents; ///< The parent IDs
      std::set<std::string> children; ///< The child IDs
      bool obsolete; ///< Flag that indicates of the term is obsolete
      std::string description; ///< Term description
      StringList synonyms; ///< List of synonyms
      StringList unparsed; ///< Unparsed lines from the definition file
      XRefType xref_type; ///< xref value-type for the CV-term
      StringList xref_binary; ///< xref binary-data-type for the CV-term (list of all allowed data value types for the current binary data array)
      std::set<std::string> units; ///< unit accession ids, defined by relationship has units

      ///Default constructor
      CVTerm();

      CVTerm(const CVTerm& rhs);

      CVTerm& operator=(const CVTerm& rhs);

      /// get mzidentml formatted string. i.e. a cvparam xml element, ref should be the name of the ControlledVocabulary (i.e. cv.name()) containing the CVTerm (e.g. PSI-MS for the psi-ms.obo - gets loaded in all cases like that??), value can be empty if not available
      std::string toXMLString(const std::string& ref, const std::string& value =std::string("")) const;

      /// get mzidentml formatted string. i.e. a cvparam xml element, ref should be the name of the ControlledVocabulary (i.e. cv.name()) containing the CVTerm (e.g. PSI-MS for the psi-ms.obo - gets loaded in all cases like that??), value can be empty if not available
      std::string toXMLString(const std::string& ref, const DataValue& value) const;

    };

    /// Constructor
    ControlledVocabulary();

    ///Destructor
    virtual ~ControlledVocabulary();

    /// Returns the CV name (set in the load method)
    const std::string& name() const;

    /// Returns the CV label (set in the load method)
    const std::string& label() const;

    /// Returns the CV version (set in the load method)
    const std::string& version() const;

    /// Returns the CV url (set in the load method)
    const std::string& url() const;

    /**
        @brief Loads the CV from an OBO file

        @param[in] name The CV name
        @param[in] filename The OBO file path

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void loadFromOBO(const std::string& name, const std::string& filename);

    /// Returns true if the term is in the CV. Returns false otherwise.
    bool exists(const std::string& id) const;

    /// Returns true if a term with the given name is in the CV. Returns false otherwise.
    bool hasTermWithName(const std::string& name) const;

    /**
        @brief Returns a term specified by ID

        @exception Exception::InvalidValue is thrown if the term is not present
    */
    const CVTerm& getTerm(const std::string& id) const;

    /**
        @brief Returns a term specified by name

        @exception Exception::InvalidValue is thrown if the term is not present
    */
    const CVTerm& getTermByName(const std::string& name, const std::string& desc = "") const;


    /// returns all the terms stored in the CV
    const std::map<std::string, CVTerm>& getTerms() const;

    /**
        @brief Writes all child terms recursively into terms

        If parent has child this method writes them recursively into the term object

        @param[out] terms Output set of child term IDs
        @param[in] parent_id The parent term ID

        @exception Exception::InvalidValue is thrown if the term is not present
    */
    void getAllChildTerms(std::set<std::string>& terms, const std::string& parent_id) const;

    /**
        @brief Writes the parent term and all descendant term IDs into @p terms

        @param[in,out] terms Output set extended with @p parent_id and all descendants
        @param[in] parent_id The parent term ID

        @exception Exception::InvalidValue is thrown if the term is not present
    */
    void addAllChildTerms(std::set<std::string>& terms, const std::string& parent_id) const;

    /**
        @brief Iterates over all children (incl. subchildren etc) of parent recursively, i.e. the whole subtree.

        @param[in] parent_id Id of parent (to be passed to getTerm(), to obtain its children).
        @param[in] lbd Function that gets the child-Ids passed. Must return bool.
                 Used for comparisons and / or to set captured variables.
                 If the lambda returns true, the iteration is exited prematurely.
                 E.g. if you have found your search, you don't need to continue searching.
                 Otherwise, if you want to go through the whole tree (e.g. to fill a vector)
                 you can just return false always to not quit early.
    */
    template <class LAMBDA>
    bool iterateAllChildren(const std::string& parent_id, LAMBDA lbd) const
    {
      for (const auto& child_id : getTerm(parent_id).children)
      {
        if (lbd(child_id) || iterateAllChildren(child_id, lbd))
          return true;
      }
      return false;
    }

    /**
        @brief Searches the existing terms for the given @p name

        @param[in] name The term name to search for

        @return const Pointer to found term. When term is not found, returns nullptr
    */
    const ControlledVocabulary::CVTerm* checkAndGetTermByName(const std::string& name) const;

    /**
        @brief Returns if @p child is a child of @p parent

        @param[in] child_id The child term ID
        @param[in] parent_id The parent term ID

        @exception Exception::InvalidValue is thrown if one of the terms is not present
    */
    bool isChildOf(const std::string& child_id, const std::string& parent_id) const;


    /**
      @brief Returns a CV for parsing/storing PSI-MS related data, e.g. mzML, or handle accessions/ids in datastructures

      The CV will be initialized on first access. Repeated access is therefor cheap.

      It consists of the following CVs:<br>
      <ul>
        <li>PSI-MS (psi-ms.obo)</li>
        <li>PATO (quality.obo)</li>
        <li>UO (unit.obo)</li>
        <li>BTO (CV/brenda.obo)</li>
        <li>GO (goslim_goa.obo)</li>
      </ul>
    */
    static const ControlledVocabulary& getPSIMSCV();

protected:
    /**
        @brief checks if a name corresponds to an id

        If the term is not known, 'true' is returned!
    */
    bool checkName_(const std::string& id, const std::string& name, bool ignore_case = true) const;

    /// Map from ID to CVTerm
    // note: unordered_map would be faster (5% for loading mzML), but order differs across platforms
    std::map<std::string, CVTerm> terms_;
    /// Map from name to id
    std::map<std::string, std::string> namesToIds_;
    /// Name set in the load method
    std::string name_;
    /// CV label
    std::string label_;
    /// CV version
    std::string version_;
    /// CV URL
    std::string url_;
  };

  ///Print the contents to a stream.
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const ControlledVocabulary& cv);


} // namespace OpenMS

