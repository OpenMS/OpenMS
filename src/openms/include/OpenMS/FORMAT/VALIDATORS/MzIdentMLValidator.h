// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>

namespace OpenMS
{
  class ControlledVocabulary;
  namespace Internal
  {

    /**
      @brief Semantic validator for mzIdentML files.

      Thin specialisation of @ref SemanticValidator that enables unit
      checking; all rule matching is performed by the base class against
      the supplied CV-mapping rules and controlled vocabularies.

      @ingroup FileIO
    */
    class OPENMS_DLLAPI MzIdentMLValidator :
      public SemanticValidator
    {
public:
      /**
        @brief Constructor.

        @param[in] mapping The mapping rules
        @param[in] cv @em All controlled vocabularies required for the mapping
      */
      MzIdentMLValidator(const CVMappings & mapping, const ControlledVocabulary & cv);

      /// Destructor
      ~MzIdentMLValidator() override;

private:

      /// Not implemented
      MzIdentMLValidator();

      /// Not implemented
      MzIdentMLValidator(const MzIdentMLValidator & rhs);

      /// Not implemented
      MzIdentMLValidator & operator=(const MzIdentMLValidator & rhs);

    };

  }   // namespace Internal

} // namespace OpenMS

