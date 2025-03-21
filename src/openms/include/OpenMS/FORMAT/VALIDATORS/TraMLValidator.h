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
      @brief Semantically validates MzXML files.
    */
    class OPENMS_DLLAPI TraMLValidator :
      public SemanticValidator
    {
public:
      /**
        @brief Constructor

                @param mapping The mapping rules
                @param cv @em All controlled vocabularies required for the mapping
            */
      TraMLValidator(const CVMappings & mapping, const ControlledVocabulary & cv);

      /// Destructor
      ~TraMLValidator() override;

private:

      /// Not implemented
      TraMLValidator();

      /// Not implemented
      TraMLValidator(const TraMLValidator & rhs);

      /// Not implemented
      TraMLValidator & operator=(const TraMLValidator & rhs);

    };

  }   // namespace Internal

} // namespace OpenMS

