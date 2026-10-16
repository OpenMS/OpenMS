// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

#include <vector>

namespace OpenMS
{
  namespace Internal
  {
    /**
    @brief Handler that is used for parsing XTandemXML data

  */
    class OPENMS_DLLAPI UnimodXMLHandler :
      public XMLHandler
    {
public:
      /// Default constructor
      UnimodXMLHandler(std::vector<ResidueModification*>& mods, const std::string& filename);

      /// Destructor
      ~UnimodXMLHandler() override;

      // Docu in base class
      void onEndElement(const char16_t* qname) override;

      // Docu in base class
      void onStartElement(const char16_t* qname, const XMLAttributes& attributes) override;

      // Docu in base class
      void onCharacters(const char16_t* chars, Size /*length*/) override;

private:

      std::string tag_;

      double avge_mass_;

      double mono_mass_;

      EmpiricalFormula diff_formula_;

      /**
          @brief Whether \<element\> children currently contribute to #diff_formula_

          Only \<delta\> (the composition of the modification itself) and \<NeutralLoss\> (the
          composition of a neutral loss) describe the chemistry of the modification. \<Ignore\>
          (reagent fragments a search must not consider), \<aa\> and \<brick\> carry
          \<element\> children as well, and folding those into the formula would leave it
          inconsistent with the recorded delta mass.
      */
      bool collect_elements_{false};

      std::vector<EmpiricalFormula> neutral_loss_diff_formula_;

      /// masses of the \<NeutralLoss\> currently open; the \<delta\> masses live in #mono_mass_ / #avge_mass_
      double neutral_loss_mono_mass_{0.0};
      double neutral_loss_avge_mass_{0.0};

      bool was_valid_peptide_modification_;
      std::vector<std::vector<EmpiricalFormula>> neutral_loss_diff_formulas_;
      std::vector<double> neutral_loss_mono_masses_;
      std::vector<double> neutral_loss_avg_masses_;
      /// one entry per specificity, parallel to #neutral_loss_diff_formulas_
      std::vector<std::vector<double>> neutral_loss_mono_masses_all_;
      std::vector<std::vector<double>> neutral_loss_avg_masses_all_;

      ResidueModification* modification_;

      std::vector<ResidueModification*>& modifications_;

      std::vector<char> sites_;

      std::vector<ResidueModification::TermSpecificity> term_specs_;
    };

  }   // namespace Internal
} // namespace OpenMS
