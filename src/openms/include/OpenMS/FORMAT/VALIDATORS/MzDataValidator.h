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
      @brief Semantically validates MzXML files.
    */
    class OPENMS_DLLAPI MzDataValidator :
      public SemanticValidator
    {
public:
      /**
        @brief Constructor

                @param mapping The mapping rules
                @param cv @em All controlled vocabularies required for the mapping
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

