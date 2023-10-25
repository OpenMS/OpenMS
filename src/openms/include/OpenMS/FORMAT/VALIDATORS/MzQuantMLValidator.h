// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#pragma once


#include <OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>

#include <map>


namespace OpenMS
{
  class ControlledVocabulary;
  namespace Internal
  {

    /**
      @brief Semantically validates MzQuantML files.
    */
    class OPENMS_DLLAPI MzQuantMLValidator :
      public SemanticValidator
    {
public:
      /**
        @brief Constructor

                @param mapping The mapping rules
                @param cv @em All controlled vocabularies required for the mapping
            */
      MzQuantMLValidator(const CVMappings & mapping, const ControlledVocabulary & cv);

      /// Destructor
      ~MzQuantMLValidator() override;

protected:
      ///CV terms which can have a value (term => value type) - see MzMLValidator impl.
      std::map<String, std::vector<CVTerm> > param_groups_;

private:

      /// Not implemented
      MzQuantMLValidator();

      /// Not implemented
      MzQuantMLValidator(const MzQuantMLValidator & rhs);

      /// Not implemented
      MzQuantMLValidator & operator=(const MzQuantMLValidator & rhs);

    };

  }   // namespace Internal

} // namespace OpenMS

