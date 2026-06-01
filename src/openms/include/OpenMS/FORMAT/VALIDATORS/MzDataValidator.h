// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>

namespace OpenMS
{
  class ControlledVocabulary;
  namespace Internal
  {

    /**
      @brief Semantic validator for mzData files.

      Specialises @ref SemanticValidator with mzData-specific checks in
      addition to the generic CV-mapping checks: for each CV term encountered,
      the recorded unit accession (if any) must be listed among the term's
      allowed units, and the term name attribute must match the name registered
      in the controlled vocabulary (whitespace is ignored and the comparison is
      case-insensitive). Mismatches are reported as errors.

      Unit checking is enabled.

      @ingroup FileIO
    */
    class OPENMS_DLLAPI MzDataValidator :
      public SemanticValidator
    {
public:
      /**
        @brief Constructor.

        @param[in] mapping The mapping rules
        @param[in] cv @em All controlled vocabularies required for the mapping
      */
      MzDataValidator(const CVMappings & mapping, const ControlledVocabulary & cv);

      /// Destructor
      ~MzDataValidator() override;

protected:

      //Docu in base class
      void handleTerm_(const String & path, const CVTerm & parsed_term) override;

private:

      /// Not implemented
      MzDataValidator();

      /// Not implemented
      MzDataValidator(const MzDataValidator & rhs);

      /// Not implemented
      MzDataValidator & operator=(const MzDataValidator & rhs);

    };

  }   // namespace Internal

} // namespace OpenMS

