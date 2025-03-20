// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/FORMAT/HANDLERS/UnimodXMLHandler.h>

using namespace std;
using namespace xercesc;

namespace OpenMS::Internal
{


    UnimodXMLHandler::UnimodXMLHandler(vector<ResidueModification*>& mods, const String& filename) :
      XMLHandler(filename, "2.0"),
      avge_mass_(0.0),
      mono_mass_(0.0),
      modification_(nullptr),
      modifications_(mods)
    {
    }

    UnimodXMLHandler::~UnimodXMLHandler()
    = default;

    void UnimodXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
    {

      tag_ = String(sm_.convert(qname));

      // new modification?
      if (tag_ == "umod:mod" || tag_ == "mod")
      {
        sites_.clear();
        modification_ = new ResidueModification();
        String title(attributeAsString_(attributes, "title"));
        modification_->setId(title);

        String full_name(attributeAsString_(attributes, "full_name"));
        // full_name.substitute("®", ""); // remove damn character (will be interpreted differently across platforms)
        // deleted this in the unimod.xml file
        modification_->setFullName(full_name);

        Int record_id(attributeAsInt_(attributes, "record_id"));
        modification_->setUniModRecordId(record_id);
        return;
      }

      // which residues are allowed?
      if (tag_ == "umod:specificity" || tag_ == "specificity")
      {
        // neutral_loss_diff_formula_ = EmpiricalFormula();
        neutral_loss_diff_formula_.clear();

        // classification of mod
        // TODO do this for all mods, do not overwrite for each specificity
        String classification(attributeAsString_(attributes, "classification"));
        modification_->setSourceClassification(classification);

        // allowed site
        String site = attributeAsString_(attributes, "site");
        //sites_.push_back(site);

        // allowed positions
        ResidueModification::TermSpecificity position = ResidueModification::ANYWHERE;
        String pos(attributeAsString_(attributes, "position"));
        if (pos == "Anywhere")
        {
          position = ResidueModification::ANYWHERE;
        }
        else if (pos == "Protein N-term")
        {
          position = ResidueModification::PROTEIN_N_TERM;
        }
        else if (pos == "Protein C-term")
        {
          position = ResidueModification::PROTEIN_C_TERM;
        }
        else if (pos == "Any C-term")
        {
          position = ResidueModification::C_TERM;
        }
        else if (pos == "Any N-term")
        {
          position = ResidueModification::N_TERM;
        }
        else
        {
          warning(LOAD, String("Don't know allowed position called: '") + pos  + "' - setting to anywhere");
        }

        was_valid_peptide_modification_ = true;
        term_specs_.push_back(position);
        if (site.size() > 1)
        {
          site = "X"; // C-term/N-term
        }
        sites_.push_back(site[0]);
        return;
      }


      if (tag_ == "umod:NeutralLoss" || tag_ == "NeutralLoss")
      {
        // mono_mass="97.976896" avge_mass="97.9952" flag="false"
        //                               composition="H(3) O(4) P">

      }

      // delta mass definitions?
      if (tag_ == "umod:delta" || tag_ == "delta")
      {
        // avge_mass="-0.9848" mono_mass="-0.984016" composition="H N O(-1)" >
        avge_mass_ = String(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("avge_mass").c_str())))).toDouble();
        mono_mass_ = String(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("mono_mass").c_str())))).toDouble();
        return;
      }

      // <umod:element symbol="H" number="1"/>
      if (tag_ == "umod:element")
      {
        String symbol = sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("symbol").c_str())));
        String num = sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("number").c_str())));
        String isotope, tmp_symbol;
        for (Size i = 0; i != symbol.size(); ++i)
        {
          if (isdigit(symbol[i]))
          {
            isotope += symbol[i];
          }
          else
          {
            tmp_symbol += symbol[i];
          }
        }

        String formula;
        if (!isotope.empty())
        {
          formula = '(' + isotope + ')' + tmp_symbol + String(num);
        }
        else
        {
          formula = tmp_symbol + num;
        }
        diff_formula_ += EmpiricalFormula(formula);
      }
    }

    void UnimodXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
      tag_ = String(sm_.convert(qname));

      // write the modifications to vector
      if (tag_ == "umod:mod" || tag_ == "mod")
      {
        modification_->setDiffAverageMass(avge_mass_);
        modification_->setDiffMonoMass(mono_mass_);
        modification_->setDiffFormula(diff_formula_);
        for (Size i = 0; i != sites_.size(); ++i)
        {
          ResidueModification* new_mod = new ResidueModification(*modification_);
          new_mod->setOrigin(sites_[i]);
          new_mod->setTermSpecificity(term_specs_[i]);
          new_mod->setNeutralLossDiffFormulas(neutral_loss_diff_formulas_[i]);
          modifications_.push_back(new_mod);
        }

        avge_mass_ = 0.0;
        mono_mass_ = 0.0;
        diff_formula_ = EmpiricalFormula();
        term_specs_.clear();
        sites_.clear();
        neutral_loss_diff_formulas_.clear();

        delete modification_;
        return;
      }

      if (tag_ == "umod:specificity" || tag_ == "specificity")
      {
        if (was_valid_peptide_modification_) // as we exclude "Protein" modifications (see above)
        {
          neutral_loss_diff_formulas_.push_back(neutral_loss_diff_formula_);
          modification_->setNeutralLossMonoMasses(neutral_loss_mono_masses_);
          modification_->setNeutralLossAverageMasses(neutral_loss_avg_masses_);
          neutral_loss_diff_formula_.clear();
          neutral_loss_mono_masses_.clear();
          neutral_loss_avg_masses_.clear();
        }
      }

      if ( (tag_ == "umod:NeutralLoss" || tag_ == "NeutralLoss") && !diff_formula_.isEmpty() )
      {
        // now diff_formula_ contains the neutral loss diff formula
        neutral_loss_diff_formula_.push_back(diff_formula_);
        neutral_loss_mono_masses_.push_back(mono_mass_);
        neutral_loss_avg_masses_.push_back(avge_mass_);
        avge_mass_ = 0.0;
        mono_mass_ = 0.0;
        diff_formula_ = EmpiricalFormula();
      }
    }

    void UnimodXMLHandler::characters(const XMLCh* const /*chars*/, const XMLSize_t /*length*/)
    {
      // nothing to do here
    }

} // namespace OpenMS    // namespace Internal
