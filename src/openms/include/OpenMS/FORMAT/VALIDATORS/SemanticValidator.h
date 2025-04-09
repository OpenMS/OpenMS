// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/DATASTRUCTURES/CVMappings.h>
#include <map>


namespace OpenMS
{
  class ControlledVocabulary;
  namespace Internal
  {

    /**
      @brief Semantically validates XML files using CVMappings and a ControlledVocabulary.

      This is the general validator.
      @n Specialized validators for specific file formats are derived from this class.
    */
    class OPENMS_DLLAPI SemanticValidator :
      protected Internal::XMLHandler,
      public Internal::XMLFile
    {
public:
      /**
        @brief Constructor

        @param mapping The mapping rules
        @param cv @em All controlled vocabularies required for the mapping
      */
      SemanticValidator(const CVMappings & mapping, const ControlledVocabulary & cv);

      /// Destructor
      ~SemanticValidator() override;

      ///Representation of a parsed CV term
      struct CVTerm
      {
        String accession;
        String name;
        String value;
        bool has_value{};
        String unit_accession;
        bool has_unit_accession{};
        String unit_name;
        bool has_unit_name{};
      };

      /**
          @brief Semantically validates an XML file.

          @param filename The file to validate
          @param errors Errors during the validation are returned in this output parameter.
          @param warnings Warnings during the validation are returned in this output parameter.

          @return @em true if the validation was successful, @em false otherwise.

          @exception Exception::FileNotFound is thrown if the file could not be opened
      */
      bool validate(const String & filename, StringList & errors, StringList & warnings);

      /// Checks if a CVTerm is allowed in a given path
      bool locateTerm(const String & path, const CVTerm & parsed_term) const;

      /// Sets the CV parameter tag name (default: 'cvParam')
      void setTag(const String & tag);

      /// Sets the name of the attribute for accessions in the CV parameter tag name (default: 'accession')
      void setAccessionAttribute(const String & accession);

      /// Sets the name of the attribute for accessions in the CV parameter tag name (default: 'name')
      void setNameAttribute(const String & name);

      /// Sets the name of the attribute for accessions in the CV parameter tag name (default: 'value')
      void setValueAttribute(const String & value);

      /**
          @brief Set if CV term value types should be check (enabled by default)

          If set to true, the xsd value types are checked, and errors are given in the cases
              - CVTerm needs value but has none
              - CVTerm has value but must not have one
              - CVTerm has value, needs value but value is of wrong type
      */
      void setCheckTermValueTypes(bool check);

      /**
          @brief Set if CV term units should be check (disabled by default)

          If set to true additional checks for CVTerms are performed:
              - CVTerm that must have a unit, but has none
              - CVTerm that has a wrong unit
      */
      void setCheckUnits(bool check);

      /// Sets the name of the unit accession attribute (default: 'unitAccession')
      void setUnitAccessionAttribute(const String & accession);

      /// Sets the name of the unit name attribute (default: 'unitName')
      void setUnitNameAttribute(const String & name);

protected:

      // Docu in base class
      void startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes) override;

      // Docu in base class
      void endElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname) override;

      // Docu in base class
      void characters(const XMLCh * const chars, const XMLSize_t /*length*/) override;

      /// Returns the current element path
      virtual String getPath_(UInt remove_from_end = 0) const;

      /// Parses the CV term accession (required), name (required) and value (optional) from the XML attributes
      virtual void getCVTerm_(const xercesc::Attributes & attributes, CVTerm & parsed_term);

      //~ forward dekl. of a inner struct/class not possible in C++ - or our Library is overtemplated
      //~ /// make a SemanticValidator::CVTerm from a ControlledVocabulary::CVTerm (without any value or unit), needed for writing only cvs at the right places in the xml (i.e. with cvmapping)
      //~ virtual void makeCVTerm_(const ControlledVocabulary::CVTerm & lc, CVTerm & parsed_term);

      /// Handling of the term
      virtual void handleTerm_(const String & path, const CVTerm & parsed_term);

      /// Reference to the mappings
      const CVMappings & mapping_;

      /// Reference to the CVs
      const ControlledVocabulary & cv_;

      /// Validation errors
      StringList errors_;

      /// Validation warnings
      StringList warnings_;

      /// List of open tags
      StringList open_tags_;

      /// Rules (location => rule)
      std::map<String, std::vector<CVMappingRule> > rules_;

      /// Fulfilled rules (location => rule ID => term ID => term count )
      /// When a tag is closed, the fulfilled rules of the current location are checked against the required rules
      /// The fulfilled rules for that location are then deleted.
      std::map<String, std::map<String, std::map<String, UInt> > > fulfilled_;


      ///@name Tag and attribute names
      //@{
      String cv_tag_;
      String accession_att_;
      String name_att_;
      String value_att_;
      String unit_accession_att_;
      String unit_name_att_;
      bool check_term_value_types_;
      bool check_units_;
      //@}

private:

      /// Not implemented
      SemanticValidator();

      /// Not implemented
      SemanticValidator(const SemanticValidator & rhs);

      /// Not implemented
      SemanticValidator & operator=(const SemanticValidator & rhs);

    };

  }   // namespace Internal

} // namespace OpenMS


