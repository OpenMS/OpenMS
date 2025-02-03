// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>
#include <map>

namespace OpenMS
{
  class ControlledVocabulary;
  namespace Internal
  {

    /**
      @brief Semantically validates MzXML files.
    */
    class OPENMS_DLLAPI MzMLValidator :
      public SemanticValidator
    {
public:
      /**
        @brief Constructor

        @param mapping The mapping rules
        @param cv @em All controlled vocabularies required for the mapping
      */
      MzMLValidator(const CVMappings & mapping, const ControlledVocabulary & cv);

      /// Destructor
      ~MzMLValidator() override;

protected:

      // Docu in base class
      void startElement(const XMLCh * const /*uri*/, const XMLCh * const /*local_name*/, const XMLCh * const qname, const xercesc::Attributes & attributes) override;

      // Docu in base class
      String getPath_(UInt remove_from_end = 0) const override;

      // Docu in base class
      void handleTerm_(const String & path, const CVTerm & parsed_term) override;

      ///CV terms which can have a value (term => value type)
      std::map<String, std::vector<CVTerm> > param_groups_;

      ///Current referenceableParamGroup identifier
      String current_id_;

      ///Binary data array name
      String binary_data_array_;
      ///Binary data array type
      String binary_data_type_;

private:

      /// Not implemented
      MzMLValidator();

      /// Not implemented
      MzMLValidator(const MzMLValidator & rhs);

      /// Not implemented
      MzMLValidator & operator=(const MzMLValidator & rhs);

    };

  }   // namespace Internal

} // namespace OpenMS

