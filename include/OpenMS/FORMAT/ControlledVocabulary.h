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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm, Andreas Bertsch, Mathias Walzer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_CONTROLLEDVOCABULARY_H
#define OPENMS_FORMAT_CONTROLLEDVOCABULARY_H

#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <set>

namespace OpenMS
{
  /**
      @brief Representation of a controlled vocabulary.

      This representation only contains the information used for parsing and validation.
      All other lines are stored in the @em unparsed member of the the CVTerm struct.

  @ingroup Format
  */
  class OPENMS_DLLAPI ControlledVocabulary
  {
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const ControlledVocabulary& cv);

public:
    /// Representation of a CV term
    struct OPENMS_DLLAPI CVTerm
    {
      /// define xsd types allowed in cv term to specify their value-type
      enum XRefType
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

      static String getXRefTypeName(XRefType type);

      String name; ///< Text name
      String id; ///< Identifier
      std::set<String> parents; ///< The parent IDs
      std::set<String> children; ///< The child IDs
      bool obsolete; ///< Flag that indicates of the term is obsolete
      String description; ///< Term description
      StringList synonyms; ///< List of synonyms
      StringList unparsed; ///< Unparsed lines from the definition file
      XRefType xref_type; ///< xref value-type for the CV-term
      StringList xref_binary; ///< xref binary-data-type for the CV-term (list of all allowed data value types for the current binary data array)
      std::set<String> units; ///< unit accession ids, defined by relationship has units

      ///Default constructor
      CVTerm();

      CVTerm(const CVTerm& rhs);

      CVTerm& operator=(const CVTerm& rhs);

      /// get mzidentml formatted string. i.e. a cvparam xml element, ref should be the name of the ControlledVocabulary (i.e. cv.name()) containing the CVTerm (e.g. PSI-MS for the psi-ms.obo - gets loaded in all cases like that??), value can be empty if not available
      String toXMLString(const String& ref, const String& value = String("")) const;

      /// get mzidentml formatted string. i.e. a cvparam xml element, ref should be the name of the ControlledVocabulary (i.e. cv.name()) containing the CVTerm (e.g. PSI-MS for the psi-ms.obo - gets loaded in all cases like that??), value can be empty if not available
      String toXMLString(const String& ref, const DataValue& value) const;

    };

    /// Constructor
    ControlledVocabulary();

    ///Destructor
    virtual ~ControlledVocabulary();

    /// Returns the CV name (set in the load method)
    const String& name() const;

    /**
        @brief Loads the CV from an OBO file

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void loadFromOBO(const String& name, const String& filename);

    /// Returns true if the term is in the CV. Returns false otherwise.
    bool exists(const String& id) const;

    /// Returns true if a term with the given name is in the CV. Returns false otherwise.
    bool hasTermWithName(const String& name) const;

    /**
        @brief Returns a term specified by ID

        @exception Exception::InvalidValue is thrown if the term is not present
    */
    const CVTerm& getTerm(const String& id) const;

    /**
        @brief Returns a term specified by name

        @exception Exception::InvalidValue is thrown if the term is not present
    */
    const CVTerm& getTermByName(const String& name, const String& desc = "") const;


    /// returns all the terms stored in the CV
    const Map<String, CVTerm>& getTerms() const;

    /**
        @brief Writes all child terms recursively into terms

        If parent has child this method writes them recursively into the term object

        @exception Exception::InvalidValue is thrown if the term is not present
    */
    void getAllChildTerms(std::set<String>& terms, const String& parent) const;

    /**
        @brief Returns if @p child is a child of @p parent

        @exception Exception::InvalidValue is thrown if one of the terms is not present
    */
    bool isChildOf(const String& child, const String& parent) const;

protected:
    /**
        @brief checks if a name corresponds to an id

        If the term is not known, 'true' is returned!
    */
    bool checkName_(const String& id, const String& name, bool ignore_case = true);

    ///Map from ID to CVTerm
    Map<String, CVTerm> terms_;
    ///Map from name to id
    Map<String, String> namesToIds_;
    ///Name set in the load method
    String name_;
  };

  ///Print the contents to a stream.
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const ControlledVocabulary& cv);


} // namespace OpenMS

#endif // OPENMS_FORMAT_CONTROLLEDVOCABULARY_H
