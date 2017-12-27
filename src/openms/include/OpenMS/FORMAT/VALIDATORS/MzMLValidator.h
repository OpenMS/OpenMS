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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_VALIDATORS_MZMLVALIDATOR_H
#define OPENMS_FORMAT_VALIDATORS_MZMLVALIDATOR_H

#include <OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>

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
      Map<String, std::vector<CVTerm> > param_groups_;

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

#endif
